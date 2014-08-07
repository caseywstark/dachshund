/*
2014, the dachshund authors.
*/

#include <cmath>
#include <cstdio>

#include <iostream>
#include <Eigen/Dense>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "timer.h"
#include "linalg.h"

#include "map.h"

//
// WFX version, x = (S + N)^{-1} d
//    A = S + N
//    b = d

double
WFX_SN::A(const int i, const int j, const bool exact)
{
    const int delta_ij = (i == j);

    const double x_i = pixels[i].x;
    const double y_i = pixels[i].y;
    const double z_i = pixels[i].z;
    const double N_ii = pixels[i].N;

    const double dx = x_i - pixels[j].x;
    const double dy = y_i - pixels[j].y;
    const double dz = z_i - pixels[j].z;

    const double x_perp_2 = (dx*dx + dy*dy);
    const double x_para_2 = dz*dz;
    double S_ij;
    if (exact) {
        S_ij = signal_covar_exact(x_perp_2, x_para_2, s_params);
    }
    else {
        S_ij = signal_covar(x_perp_2, x_para_2, s_params);
    }

    const double A_ij = S_ij + N_ii * delta_ij;
    return A_ij;
}

// Just a stupid Jacobi for now
void
WFX_SN::preconditioner_solve(const double * const b, double * const x)
{
    const int n = num_pixels;
    const double var_s = s_params->var_s;

    for (int i = 0; i < n; ++i) {
        // no coord lookup for diagonal. Apc_ii = var_s + N_i
        double Apc_ii = var_s + pixels[i].N;
        x[i] = b[i] / Apc_ii;
    }
}

void
WFX_SN::A_product(const double * const a, double * const b)
{
    const int n = num_pixels;

    int i, j;
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        b[i] = 0.0;
        const double x_i = pixels[i].x;
        const double y_i = pixels[i].y;
        const double z_i = pixels[i].z;
        const double N_ii = pixels[i].N;
        for (j = 0; j < n; ++j) {
            const int delta_ij = (i == j);
            const double dx = x_i - pixels[j].x;
            const double dy = y_i - pixels[j].y;
            const double dz = z_i - pixels[j].z;
            const double x_perp_2 = (dx*dx + dy*dy);
            const double x_para_2 = dz*dz;
            const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
            const double A_ij = S_ij + N_ii * delta_ij;
            b[i] += A_ij * a[j];
        }
    }
}

void
WFX_SN::A_residual(const double * const x, double * const r)
{
    const int n = num_pixels;

    int i, j;
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        double sum_Ax_i = 0.0;
        const double x_i = pixels[i].x;
        const double y_i = pixels[i].y;
        const double z_i = pixels[i].z;
        const double N_ii = pixels[i].N;
        const double d_i = pixels[i].d;

        for (j = 0; j < n; ++j) {
            const int delta_ij = (i == j);
            const double dx = x_i - pixels[j].x;
            const double dy = y_i - pixels[j].y;
            const double dz = z_i - pixels[j].z;
            const double x_perp_2 = (dx*dx + dy*dy);
            const double x_para_2 = dz*dz;
            const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
            const double A_ij = S_ij + N_ii * delta_ij;
            sum_Ax_i += A_ij * x[j];
        }
        // r_i = b_i - sum_j A_ij x_j
        r[i] = d_i - sum_Ax_i;
    }
}

double
WFX_SN::norm_b()
{
    double sum_b2 = 0.0;
    for (int i = 0; i < num_pixels; ++i) {
        const double bi = pixels[i].d;
        sum_b2 += bi * bi;
    }
    return sqrt(sum_b2);
}

void
WFX_SN::pcg_step(double * const d, double * const q, double * const r,
    double * const s, double * const x, const bool full_residual,
    double * delta)
{
    const int n = num_pixels;

    // q = A d
    A_product(d, q);
    // alpha = delta / (d^T q)
    double alpha = *delta / ds_vector_dot(n, d, q);
    // Solution update. x = x + alpha d
    ds_vector_add_cv(n, alpha, d, x);

    // Update residual.
    if (full_residual) {
        // Full update, r = b - A x
        A_residual(x, r);
    }
    else {
        // Approximate update, r = r - alpha q
        ds_vector_sub_cv(n, alpha, q, r);
    }

    // reapply preconditioner, s = M^-1 r
    preconditioner_solve(r, s);

    // save current delta.
    double delta_old = *delta;

    // Update delta, delta_new = r^T s
    *delta = ds_vector_dot(n, r, s);

    // Update d, d = s + beta d
    double beta = *delta / delta_old;
    for (int i = 0; i < n; ++i) {
        d[i] = s[i] + beta * d[i];
    }
}

PCGResult
WFX_SN::solve_pcg(double * const x, const PCGParams * const pcg_params,
    const bool verbose)
{
    const int n = num_pixels;
    const int max_iter = pcg_params->max_iter;
    const double tol = pcg_params->tol;
    const int step_r = pcg_params->step_r;

    // return val.
    PCGResult res;

    // Work arrays
    double *d = new double[n];
    double *r = new double[n];
    double *s = new double[n];
    double *q = new double[n];

    // Compute norm b for stop condition.
    double tol_norm_b = tol * norm_b();

    // Setup the residual, r = b - A x
    A_residual(x, r);

    // Preconditioning, d = M^-1 r
    preconditioner_solve(r, d);

    // delta = r^T d
    double delta = ds_vector_dot(n, r, d);

    // residual norm
    double sum_r2 = 0.0;
    for (int i = 0; i < n; ++i) { sum_r2 += r[i] * r[i]; }
    double norm_r = sqrt(sum_r2);
    res.residual_norms.push_back(norm_r);

    if (verbose) {
        printf("[PCG] Solving %i x %i problem.\n", n, n);
        printf("    Goal |r| < (tol) |b| = %e,  delta = %e\n", tol_norm_b, delta);
    }

    int iter;
    for (iter = 1; iter < max_iter; ++iter) {
        bool full_residual = iter % step_r == 0;
        pcg_step(d, q, r, s, x, full_residual, &delta);

        // compute residual norm
        sum_r2 = 0.0;
        for (int i = 0; i < n; ++i) { sum_r2 += r[i] * r[i]; }
        norm_r = sqrt(sum_r2);
        res.residual_norms.push_back(norm_r);

        // Check stop condition.
        if (norm_r <= tol_norm_b) {
            if (verbose) {
                printf("    iter %i reached |r| = %e < tol |b| = %e\n", iter, norm_r, tol_norm_b);
            }
            break;
        }

        // Output progress.
        if (verbose) {
            printf("    iter %i, |r| %e, delta %e\n", iter, norm_r, delta);
        }
    }

    // don't forget to free the work arrays.
    delete [] d;
    delete [] r;
    delete [] q;
    delete [] s;

    // max iter warning
    if (iter == max_iter) {
        printf("[WARNING] PCG reached max iteration before stop condition.\n");
    }

    // Update res before returning
    res.num_iters = iter;
    return res;
}

void
WFX_SN::solve_cf(double * const x, const bool verbose, const bool exact)
{
    using namespace Eigen;
    const int n = num_pixels;

    if (verbose) {
        puts("Filling A matrix.");
    }

    // init A = (S + N)
    MatrixXd sn(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            sn(i, j) = A(i, j, exact);
        }
    }

    if (verbose) {
        puts("Filling b vector.");
    }

    // init b = d
    VectorXd b(n);
    for (int i = 0; i < n; ++i) {
        b(i) = pixels[i].d;
    }

    if (verbose) {
        puts("Calling eigen solve.");
    }

    // solve
    VectorXd vec_x(n);
    vec_x = sn.ldlt().solve(b);

    if (verbose) {
        puts("Copying back.");
    }

    // write back to our format.
    for (int i = 0; i < n; ++i) {
        x[i] = vec_x(i);
    }
}

//
// WFX version x = w (w S w + I)^{-1} w d
//   A_ij = w_i w_j S_ij + delta_ij
//   b_i = w_i d_i
//

double
WFX_W::A(const int i, const int j, const bool exact)
{
    const int delta_ij = (i == j);

    const double x_i = pixels[i].x;
    const double y_i = pixels[i].y;
    const double z_i = pixels[i].z;
    const double w_i = pixels[i].w;

    const double dx = x_i - pixels[j].x;
    const double dy = y_i - pixels[j].y;
    const double dz = z_i - pixels[j].z;
    const double w_j = pixels[j].w;

    const double x_perp_2 = (dx*dx + dy*dy);
    const double x_para_2 = dz*dz;

    double S_ij;
    if (exact) {
        S_ij = signal_covar_exact(x_perp_2, x_para_2, s_params);
    }
    else {
        S_ij = signal_covar(x_perp_2, x_para_2, s_params);
    }

    const double A_ij = w_i * S_ij * w_j + delta_ij;
    return A_ij;
}

void
WFX_W::preconditioner_solve(const double * const b, double * const x)
{
    const int n = num_pixels;
    const double var_s = s_params->var_s;

    for (int i = 0; i < n; ++i) {
        // no coord lookup for diagonal. Apc_ii = var_s + N_i
        double wi = pixels[i].w;
        double Apc_ii = wi * wi * var_s + 1.0;
        x[i] = b[i] / Apc_ii;
    }
}

void
WFX_W::A_product(const double * const a, double * const b)
{
    const int n = num_pixels;

    int i, j;
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        b[i] = 0.0;
        const double x_i = pixels[i].x;
        const double y_i = pixels[i].y;
        const double z_i = pixels[i].z;
        const double w_i = pixels[i].w;
        for (j = 0; j < n; ++j) {
            const int delta_ij = (i == j);
            const double w_j = pixels[j].w;
            const double dx = x_i - pixels[j].x;
            const double dy = y_i - pixels[j].y;
            const double dz = z_i - pixels[j].z;
            const double x_perp_2 = (dx*dx + dy*dy);
            const double x_para_2 = dz*dz;
            const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
            const double A_ij = w_i * S_ij * w_j + delta_ij;
            b[i] += A_ij * a[j];
        }
    }
}

void
WFX_W::A_residual(const double * const x, double * const r)
{
    const int n = num_pixels;

    int i, j;
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        double sum_Ax_i = 0.0;
        const double x_i = pixels[i].x;
        const double y_i = pixels[i].y;
        const double z_i = pixels[i].z;
        const double w_i = pixels[i].w;

        for (j = 0; j < n; ++j) {
            const int delta_ij = (i == j);
            const double w_j = pixels[j].w;
            const double dx = x_i - pixels[j].x;
            const double dy = y_i - pixels[j].y;
            const double dz = z_i - pixels[j].z;
            const double x_perp_2 = (dx*dx + dy*dy);
            const double x_para_2 = dz*dz;
            const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
            const double A_ij = w_i * S_ij * w_j + delta_ij;
            sum_Ax_i += A_ij * x[j];
        }

        // r_i = b_i - sum_j A_ij x_j
        const double wd_i = pixels[i].wd;
        r[i] = wd_i - sum_Ax_i;
    }
}

double
WFX_W::norm_b()
{
    double sum_b2 = 0.0;
    for (int i = 0; i < num_pixels; ++i) {
        const double bi = pixels[i].wd;
        sum_b2 += bi * bi;
    }
    return sqrt(sum_b2);
}

void
WFX_W::pcg_step(double * const d, double * const q, double * const r,
    double * const s, double * const x, const bool full_residual,
    double * delta)
{
    const int n = num_pixels;

    // q = A d
    A_product(d, q);
    // alpha = delta / (d^T q)
    double alpha = *delta / ds_vector_dot(n, d, q);
    // Solution update. x = x + alpha d
    ds_vector_add_cv(n, alpha, d, x);

    // Update residual.
    if (full_residual) {
        // Full update, r = b - A x
        A_residual(x, r);
    }
    else {
        // Approximate update, r = r - alpha q
        ds_vector_sub_cv(n, alpha, q, r);
    }

    // reapply preconditioner, s = M^-1 r
    preconditioner_solve(r, s);

    // save current delta.
    double delta_old = *delta;

    // Update delta, delta_new = r^T s
    *delta = ds_vector_dot(n, r, s);

    // Update d, d = s + beta d
    double beta = *delta / delta_old;
    for (int i = 0; i < n; ++i) {
        d[i] = s[i] + beta * d[i];
    }
}


PCGResult
WFX_W::solve_pcg(double * const x, const PCGParams * const pcg_params,
    const bool verbose)
{
    const int n = num_pixels;
    const int max_iter = pcg_params->max_iter;
    const double tol = pcg_params->tol;
    const int step_r = pcg_params->step_r;

    // return val.
    PCGResult res;

    // Work arrays
    double *d = new double[n];
    double *r = new double[n];
    double *s = new double[n];
    double *q = new double[n];

    // Compute norm b for stop condition.
    double tol_norm_b = tol * norm_b();

    // Setup the residual, r = b - A x
    A_residual(x, r);

    // Preconditioning, d = M^-1 r
    preconditioner_solve(r, d);

    // delta = r^T d
    double delta = ds_vector_dot(n, r, d);

    // residual norm
    double norm_r = sqrt(ds_vector_dot(n, r, r));
    res.residual_norms.push_back(norm_r);

    if (verbose) {
        printf("[PCG] Solving %i x %i problem.\n", n, n);
        printf("    Goal |r| < (tol) |b| = %e,  delta = %e\n", tol_norm_b, delta);
    }

    int iter;
    for (iter = 1; iter < max_iter; ++iter) {
        bool full_residual = iter % step_r == 0;
        pcg_step(d, q, r, s, x, full_residual, &delta);

        // compute residual norm
        norm_r = sqrt(ds_vector_dot(n, r, r));
        res.residual_norms.push_back(norm_r);

        // Check stop condition.
        if (norm_r <= tol_norm_b) {
            if (verbose) {
                printf("    iter %i reached |r| = %e < tol |b| = %e\n", iter, norm_r, tol_norm_b);
            }
            break;
        }

        // Output progress.
        if (verbose) {
            printf("    iter %i, |r| %e, delta %e\n", iter, norm_r, delta);
        }
    }

    // don't forget to free the work arrays.
    delete [] d;
    delete [] r;
    delete [] q;
    delete [] s;

    // max iter warning
    if (iter == max_iter) {
        printf("[WARNING] PCG reached max iteration before stop condition.\n");
    }

    // Don't forget the last multiply!
    for (int i = 0; i < n; ++i) {
        x[i] *= pixels[i].w;
    }

    // Update res before returning
    res.num_iters = iter;
    return res;
}

void
WFX_W::solve_cf(double * const x, const bool verbose, const bool exact)
{
    using namespace Eigen;
    const int n = num_pixels;

    if (verbose) {
        puts("Filling A matrix.");
    }

    // init A = (S + N)
    MatrixXd sn(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            sn(i, j) = A(i, j, exact);
        }
    }

    if (verbose) {
        puts("Filling b vector.");
    }

    // init b = d
    VectorXd b(n);
    for (int i = 0; i < n; ++i) {
        b(i) = pixels[i].wd;
    }

    if (verbose) {
        puts("Calling eigen solve.");
    }

    // solve
    VectorXd vec_x(n);
    vec_x = sn.ldlt().solve(b);

    if (verbose) {
        puts("Copying back.");
    }

    // write back to our format
    // Don't forget the last multiply!
    for (int i = 0; i < n; ++i) {
        x[i] = pixels[i].w * vec_x(i);
    }
}


//
// Map
//

void
smp_product(const int num_pixel_points, const Point * const pixel_coords,
        const int num_map_points, const Point * const map_coords,
        const SignalCovarParams * const s_params,
        const double * const x, double * const m)
{
    int i, j;
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < num_map_points; ++i) {
        m[i] = 0.0;
        const double x_i = map_coords[i].x;
        const double y_i = map_coords[i].y;
        const double z_i = map_coords[i].z;
        for (j = 0; j < num_pixel_points; ++j) {
            const double dx = x_i - pixel_coords[j].x;
            const double dy = y_i - pixel_coords[j].y;
            const double dz = z_i - pixel_coords[j].z;
            const double x_perp_2 = (dx*dx + dy*dy);
            const double x_para_2 = dz*dz;
            const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
            m[i] += S_ij * x[j];
        }
    }
}

/*
void
Map::solve(double * const m, const bool verbose)
{
    // Start timing.
    const int t_icat = g_timer.add_category("map");

    // Array for the pcg solution.
    double * const x = new double[num_pixels];
    // Zero out initial guess.
    for (int i = 0; i < num_pixels; ++i) { x[i] = 0.0; }

    // x = (W S_{pp} + I)^{-1} W d
    if (verbose) {
        printf("Entering PCG with num_pixels = %i\n", num_pixels);
    }
    int it0 = g_timer.record(t_icat);
    PCGResult res = x_solve(x, verbose);
    int it1 = g_timer.record(t_icat);
    if (verbose) {
        printf("[TIME] %g s for PCG.\n", g_timer.sec_interval(t_icat, it0, it1));
    }

    // m = S^{mp} x
    if (verbose) {
        printf("Entering map multiply with num_points = %i, num_points x num_pixels = %.2e.\n", num_points, (double)num_points * num_pixels);
    }
    int it2 = g_timer.record(t_icat);
    Smp_x_mult(x, m);
    int it3 = g_timer.record(t_icat);
    if (verbose) {
        printf("[TIME] %g s for map multiply.\n", g_timer.sec_interval(t_icat, it2, it3));
    }

    delete [] x;
}


void
WFX::write_A(const std::string output_path)
{
    const int n = num_pixels;

    double *a = new double[n];
    FILE *output_file = fopen(output_path.c_str(), "w");

    int i, j;
    for (i = 0; i < n; ++i) {
        const double x_i = pixels[i].x;
        const double y_i = pixels[i].y;
        const double z_i = pixels[i].z;
        const double N_ii = pixels[i].N;
        for (j = 0; j < n; ++j) {
            const int delta_ij = (i == j);
            const double dx = x_i - pixels[j].x;
            const double dy = y_i - pixels[j].y;
            const double dz = z_i - pixels[j].z;
            const double x_perp_2 = (dx*dx + dy*dy);
            const double x_para_2 = dz*dz;
            const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
            const double A_ij = S_ij + N_ii * delta_ij;
            a[j] = A_ij;
        }

        fwrite(a, sizeof(double), n, output_file);
    }

    fclose(output_file);

    delete [] a;
}


// ===
// Section for S^{mp} W S^{pm}
// ===

void
ds_smp_w_spm(double * const m)
{
    const int map_n = p.map_n;
    const int pix_n = p.num_pixels;
    const double var_s = p.corr_var_s;
    const double l_perp_2 = p.corr_l_perp_2;
    const double l_para_2 = p.corr_l_para_2;
    const int gt_n = gt.n;
    const double gt_dx = gt.dx;
    const double * const gt_table = gt.table;
    const DSPoint * const l_mc = map_coords;
    const DSPixel * const l_pixels = pixels;

    int i, j;

#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < map_n; ++i) {
        // Reset value at this map point.
        m[i] = 0.0;
        // Grab map coords.
        double xi = l_mc[i].x;
        double yi = l_mc[i].y;
        double zi = l_mc[i].z;

        for (j = 0; j < pix_n; ++j) {
            double dx = xi - l_pixels[j].x;
            double dy = yi - l_pixels[j].y;
            double dz = zi - l_pixels[j].z;
            double wj = l_pixels[j].w;
            double x_perp_2 = (dx*dx + dy*dy) / l_perp_2;
            double x_para_2 = dz*dz / l_para_2;
            double g_perp = ds_linterp(x_perp_2, gt_n, gt_dx, gt_table);
            double g_para = ds_linterp(x_para_2, gt_n, gt_dx, gt_table);
            double smp = var_s * g_perp * g_para;
            m[i] += smp * wj * smp;
        }
    }
}

// ===
// Section for (W S^{pp} + I)^{-1} S^{pm}
// ===

void
ds_wsppi_spm_residual(const int k, const double * const x,
    double * const r)
{
    const double xk = map_coords[k].x;
    const double yk = map_coords[k].y;
    const double zk = map_coords[k].z;
    const double var_s = p.corr_var_s;
    const double l_para_2 = p.corr_l_para_2;
    const double l_perp_2 = p.corr_l_perp_2;
    const int gt_n = gt.n;
    const double gt_dx = gt.dx;
    const double * const gt_table = gt.table;
    const DSPixel * const l_pixels = pixels;

    const int n = p.num_pixels;
    int i, j;

#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        double Ax_i = 0.0;
        double xi = l_pixels[i].x;
        double yi = l_pixels[i].y;
        double zi = l_pixels[i].z;
        double one_minus_wi = 1.0 - l_pixels[i].w;
        for (j = 0; j < n; ++j) {
            double dx = xi - l_pixels[j].x;
            double dy = yi - l_pixels[j].y;
            double dz = zi - l_pixels[j].z;
            double x_perp_2 = (dx*dx + dy*dy) / l_perp_2;
            double x_para_2 = dz*dz / l_para_2;
            double g_perp = ds_linterp(x_perp_2, gt_n, gt_dx, gt_table);
            double g_para = ds_linterp(x_para_2, gt_n, gt_dx, gt_table);
            double wsppi_ij = var_s * g_perp * g_para;
            int delta_ij = (i == j);
            wsppi_ij = wsppi_ij + (1.0 - one_minus_wi * wsppi_ij) * delta_ij;
            Ax_i += wsppi_ij * x[j];
        }
        double dx = xi - xk;
        double dy = yi - yk;
        double dz = zi - zk;
        double x_perp_2 = (dx*dx + dy*dy) / l_perp_2;
        double x_para_2 = dz*dz / l_para_2;
        double g_perp = ds_linterp(x_perp_2, gt_n, gt_dx, gt_table);
        double g_para = ds_linterp(x_para_2, gt_n, gt_dx, gt_table);
        double spm_ik = var_s * g_perp * g_para;
        r[i] = spm_ik - Ax_i;
    }
}

void
ds_wsppi_spm_q(const double * const d, double * const q)
{
    const double var_s = p.corr_var_s;
    const double l_para_2 = p.corr_l_para_2;
    const double l_perp_2 = p.corr_l_perp_2;
    const int gt_n = gt.n;
    const double gt_dx = gt.dx;
    const double * const gt_table = gt.table;
    const DSPixel * const l_pixels = pixels;

    const int n = p.num_pixels;
    int i, j;

#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        q[i] = 0.0;
        double xi = l_pixels[i].x;
        double yi = l_pixels[i].y;
        double zi = l_pixels[i].z;
        double one_minus_wi = 1.0 - l_pixels[i].w;
        for (j = 0; j < n; ++j) {
            double dx = xi - l_pixels[j].x;
            double dy = yi - l_pixels[j].y;
            double dz = zi - l_pixels[j].z;
            double x_perp_2 = (dx*dx + dy*dy) / l_perp_2;
            double x_para_2 = dz*dz / l_para_2;
            double g_perp = ds_linterp(x_perp_2, gt_n, gt_dx, gt_table);
            double g_para = ds_linterp(x_para_2, gt_n, gt_dx, gt_table);
            double wsppi_ij = var_s * g_perp * g_para;
            int delta_ij = (i == j);
            wsppi_ij = wsppi_ij + (1.0 - one_minus_wi * wsppi_ij) * delta_ij;
            q[i] += wsppi_ij * d[j];
        }
    }
}

void
pcg_wsppi_spm(const int k, double * const x, const double tol_norm_b,
    double * const d, double * const r, double * const s, double * const q)
{
    int i;
    const int n = p.num_pixels;
    const int max_iter = p.pcg_max_iter;
    const int pcg_step_r = p.pcg_step_r;
    const double corr_var_s = p.corr_var_s;

    // Setup the residual.
    // r = b - A x
    ds_wsppi_spm_residual(k, x, r);

    // The preconditioning.
    // For now, we will just try a Jacobian M (A diag).
    // s = M^-1 r
    for (i = 0; i < n; ++i) {
        double wsppi_i = pixels[i].w * corr_var_s + 1.0;
        d[i] = r[i] / wsppi_i;
    }

    // delta = r^T d
    double delta_new = 0.0;
    for (i = 0; i < n; ++i) {
        delta_new += r[i] * d[i];
    }

    //printf("    delta_0 %e, tol %e\n", delta_new, tol_norm_b);

    for (int iter = 1; iter < max_iter; ++iter) {
        // q = A d
        ds_wsppi_spm_q(d, q);

        // alpha = delta_new / (d^T q)
        double denom = 0.0;
        for (i = 0; i < n; ++i) {
            denom += d[i] * q[i];
        }
        double alpha = delta_new / denom;

        // x = x + alpha d
        for (i = 0; i < n; ++i) {
            x[i] += alpha * d[i];
        }

        // Update residual.
        if (iter % pcg_step_r == 0) {
            // Full update.
            // r = b - A x
            ds_wsppi_spm_residual(k, x, r);
        }
        else {
            // Approximate update.
            // r = r - alpha q
            for (i = 0; i < n; ++i) {
                r[i] -= alpha * q[i];
            }
        }

        // reapply preconditioner.
        for (i = 0; i < n; ++i) {
            double wsppi_ii = pixels[i].w * corr_var_s + 1.0;
            s[i] = r[i] / wsppi_ii;
        }

        // save current delta.
        double delta_old = delta_new;

        // Update delta.
        // delta_new = r^T s
        delta_new = 0.0;
        for (i = 0; i < n; ++i) {
            delta_new += r[i] * s[i];
        }

        // Update d.
        double beta = delta_new / delta_old;
        for (i = 0; i < n; ++i) {
            d[i] = s[i] + beta * d[i];
        }

        // Check stop condition.
        if (delta_new <= tol_norm_b) {
            //printf("    reached tolerance %e, on iter %i, delta %e\n", tol_norm_b, iter, delta_new);
            break;
        }

        // Output progress.
        //printf("    iter %i, delta %e\n", iter, delta_new);
    }
}

void
ds_map_covar_diag(double * const M_diag)
{
    int i, j;

    // constant local stuff.
    const int pix_n = p.num_pixels;
    const int map_n = p.map_n;
    const double tol = p.pcg_tol;
    const double var_s = p.corr_var_s;
    const double l_perp_2 = p.corr_l_perp_2;
    const double l_para_2 = p.corr_l_para_2;
    const int gt_n = gt.n;
    const double gt_dx = gt.dx;
    const double * const gt_table = gt.table;
    const DSPixel * const l_pixels = pixels;

    // work arrays
    double * const x = new double[pix_n];
    double * const d = new double[pix_n];
    double * const r = new double[pix_n];
    double * const s = new double[pix_n];
    double * const q = new double[pix_n];

    // Iterate over map points.
    for (i = 0; i < map_n; ++i) {
        // Print progress because this part is sloooow.
        printf("cell %i -- %f.\n", i, (double)i/map_n);

        // Reset M val.
        M_diag[i] = 0.0;

        // Compute norm b for stop condition.
        double norm_b = 0.0;
        const double xi = map_coords[i].x;
        const double yi = map_coords[i].y;
        const double zi = map_coords[i].z;

        for (j = 0; j < pix_n; ++j) {
            double dx = xi - l_pixels[j].x;
            double dy = yi - l_pixels[j].y;
            double dz = zi - l_pixels[j].z;
            double x_perp_2 = (dx*dx + dy*dy) / l_perp_2;
            double x_para_2 = dz*dz / l_para_2;
            double g_perp = ds_linterp(x_perp_2, gt_n, gt_dx, gt_table);
            double g_para = ds_linterp(x_para_2, gt_n, gt_dx, gt_table);
            double spm_ik = var_s * g_perp * g_para;
            norm_b += spm_ik * spm_ik;
        }

        norm_b = sqrt(norm_b);
        double tol_norm_b = tol * norm_b;

        // Annoying edge case...
        // Sometimes [Sum_i Spm_ik^2]^{1/2} = 0, meaning that the map point has
        // no overlap with the pixels. In this case, M_ii = 0 anyway, so don't
        // continue.
        if (norm_b > 0.0) {
            // PCG for x_j = (W S^{pp} + I)^{-1}_jk S^{pm}_ki
            pcg_wsppi_spm(i, x, tol_norm_b, d, r, q, s);

            // Next step M_{ii} = S^{mp}_{ij} x_j
            for (j = 0; j < pix_n; ++j) {
                // Separation vector components.
                double dx = xi - l_pixels[j].x;
                double dy = yi - l_pixels[j].y;
                double dz = zi - l_pixels[j].z;
                double x_perp_2 = (dx*dx + dy*dy) / l_perp_2;
                double x_para_2 = dz*dz / l_para_2;
                double g_perp = ds_linterp(x_perp_2, gt_n, gt_dx, gt_table);
                double g_para = ds_linterp(x_para_2, gt_n, gt_dx, gt_table);
                double smp_ij = var_s * g_perp * g_para;
                M_diag[i] += smp_ij * x[j];
            }
        }
    }

    // don't forget to free the work arrays.
    delete [] x;
    delete [] d;
    delete [] r;
    delete [] q;
    delete [] s;
}
*/
