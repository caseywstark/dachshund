/*
2014, the dachshund authors.
*/

#include <cmath>
#include <cstdio>
#include <omp.h>

#include "params.h"
#include "timer.h"

#include "pcg.h"

struct ds_point_s {
    double x;
    double y;
    double z;
};

typedef struct ds_point_s DSPoint;

struct ds_pixel_s {
    double x;
    double y;
    double z;
    double w;
};

typedef struct ds_pixel_s DSPixel;

struct ds_func_table_s {
    int n;
    double dx;
    double *table;
};

typedef struct ds_func_table_s DSFuncTable;

// global static variables
DSFuncTable gt;
DSPixel *pixels;
DSPoint *map_coords;

void
ds_gt_init(const int n, const double dx)
{
    gt.n = n;
    gt.dx = dx;
    gt.table = new double[n];
    for (int i = 0; i < n; ++i) {
        double x = i * dx;
        gt.table[i] = exp(-x);
    }
}

void
ds_gt_free()
{
    delete gt.table;
}

void
ds_pixels_init(const double * const pixel_x, const double * const pixel_y,
    const double * const pixel_z, const double * const pixel_w)
{
    const int pix_n = p.num_pixels;
    pixels = new DSPixel[pix_n];

    for (int i = 0; i < pix_n; ++i) {
        pixels[i].x = pixel_x[i];
        pixels[i].y = pixel_y[i];
        pixels[i].z = pixel_z[i];
        pixels[i].w = pixel_w[i];
    }
}

void
ds_pixels_free()
{
    delete pixels;
}

void
ds_map_coords_init(const int map_nx, const int map_ny, const int map_nz,
    const double map_dx, const double map_dy, const double map_dz)
{
    int map_n = map_nx * map_ny * map_nz;
    map_coords = new DSPoint[map_n];

    for (int ix = 0; ix < map_nx; ++ix) {
        for (int iy = 0; iy < map_ny; ++iy) {
            for (int iz = 0; iz < map_nz; ++iz) {
                int i = (ix * map_ny + iy) * map_nz + iz;
                map_coords[i].x = map_dx * (ix + 0.5);
                map_coords[i].y = map_dy * (iy + 0.5);
                map_coords[i].z = map_dz * (iz + 0.5);
            }
        }
    }
}

void
ds_map_coords_free()
{
    delete map_coords;
}

void
ds_lum_init(const double * const pixel_x, const double * const pixel_y,
    const double * const pixel_z,  const double * const pixel_w)
{
    ds_pixels_init(pixel_x, pixel_y, pixel_z, pixel_w);
    ds_map_coords_init(p.map_nx, p.map_ny, p.map_nz, p.map_dx, p.map_dy,
        p.map_dz);

    int gt_n = 1000;
    double gt_dx = 10.0 / gt_n;
    ds_gt_init(gt_n, gt_dx);
}

void
ds_lum_free()
{
    ds_map_coords_free();
    ds_pixels_free();
    ds_gt_free();
}

/*
Linear interpolation, given x and table data.

We assume x >= 0. If x is larger than the last tabulated point, we return 0.

Parameters
x : the eval point
nx : the number of samples
dx: the sample point separation
f_table : array containing the function value at sample points.
*/
static inline double
ds_linterp(const double x, const int nx, const double dx,
    const double * const f_table)
{
    // in index space.
    double xx = x / dx;
    // decimate to get low edge.
    int ix = xx;
    // Check if we are outside table range.
    if (ix > nx - 2) {
        return 0.0;
    }
    // the linterp expression.
    return f_table[ix] * (ix + 1 - xx) + f_table[ix+1] * (xx - ix);
}

// ===
// Section for (W S^{pp} + I)^{-1} b
// ===

void
ds_wsppi_residual(const double * const x, const double * const b,
    double * const r)
{
    // Make local copies of global vars for access performance.
    const int n = p.num_pixels;
    const double var_s = p.corr_var_s;
    const double l_perp_2 = p.corr_l_perp_2;
    const double l_para_2 = p.corr_l_para_2;
    const int gt_n = gt.n;
    const double gt_dx = gt.dx;
    const double * const gt_table = gt.table;
    const DSPixel * const l_pixels = pixels;

    int i, j;
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        double Ax_i = 0.0;

        double xi = l_pixels[i].x;
        double yi = l_pixels[i].y;
        double zi = l_pixels[i].z;
        double wii = l_pixels[i].w;

        for (j = 0; j < n; ++j) {
            // Separation vector components.
            double dx = xi - l_pixels[j].x;
            double dy = yi - l_pixels[j].y;
            double dz = zi - l_pixels[j].z;

            // The gaussian terms.
            double x_perp_2 = (dx*dx + dy*dy) / l_perp_2;
            double x_para_2 = dz*dz / l_para_2;

            // S = sigma^2 exp(...) exp(...)
            double g_perp = ds_linterp(x_perp_2, gt_n, gt_dx, gt_table);
            double g_para = ds_linterp(x_para_2, gt_n, gt_dx, gt_table);
            double wsppi_ij = wii * var_s * g_perp * g_para;

            // nice trick to avoid branching.
            int delta_ij = (i == j);
            wsppi_ij = wsppi_ij + delta_ij;

            // Update sum.
            Ax_i += wsppi_ij * x[j];
        }
        r[i] = b[i] - Ax_i;
    }
}

void
ds_wsppi_q(const double * const d, double * const q)
{
    // Make local copies of global vars for access performance.
    const int n = p.num_pixels;
    const double var_s = p.corr_var_s;
    const double l_perp_2 = p.corr_l_perp_2;
    const double l_para_2 = p.corr_l_para_2;
    const int gt_n = gt.n;
    const double gt_dx = gt.dx;
    const double * const gt_table = gt.table;
    const DSPixel * const l_pixels = pixels;

    int i, j;
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        q[i] = 0.0;
        double xi = l_pixels[i].x;
        double yi = l_pixels[i].y;
        double zi = l_pixels[i].z;
        double wii = l_pixels[i].w;
        for (j = 0; j < n; ++j) {
            double dx = xi - l_pixels[j].x;
            double dy = yi - l_pixels[j].y;
            double dz = zi - l_pixels[j].z;
            double x_perp_2 = (dx*dx + dy*dy) / l_perp_2;
            double x_para_2 = dz*dz / l_para_2;
            double g_perp = ds_linterp(x_perp_2, gt_n, gt_dx, gt_table);
            double g_para = ds_linterp(x_para_2, gt_n, gt_dx, gt_table);
            double wsppi_ij = wii * var_s * g_perp * g_para;
            int delta_ij = (i == j);
            wsppi_ij = wsppi_ij + delta_ij;
            q[i] += wsppi_ij * d[j];
        }
    }
}

void
pcg_wsppi_b(const double * const b, double * const x)
{
    int i;
    const int n = p.num_pixels;
    const int max_iter = p.pcg_max_iter;
    const double tol = p.pcg_tol;
    const int pcg_step_r = p.pcg_step_r;
    const double corr_var_s = p.corr_var_s;

    // Work arrays
    double *d = new double[n];
    double *r = new double[n];
    double *s = new double[n];
    double *q = new double[n];

    // Compute norm b for stop condition.
    double norm_b = 0.0;
    for (i = 0; i < n; ++i) {
        norm_b += b[i]*b[i];
    }
    norm_b = sqrt(norm_b);
    double tol_norm_b = tol * norm_b;

    // Setup the residual.
    // r = b - A x
    ds_wsppi_residual(x, b, r);

    // The preconditioning.
    // For now, we will just try a Jacobian M (A diag).
    // s = M^-1 r
    for (i = 0; i < n; ++i) {
        // no coord lookup for diagonal piece...
        // just w_i var + 1
        double wsppi_ii = pixels[i].w * corr_var_s + 1.0;
        d[i] = r[i] / wsppi_ii;
    }

    // delta = r^T d
    double delta_new = 0.0;
    for (i = 0; i < n; ++i) {
        delta_new += r[i] * d[i];
    }

    printf("    delta_0 = %e, tol |b| = %e\n", delta_new, tol_norm_b);

    for (int iter = 1; iter < max_iter; ++iter) {
        // q = A d
        ds_wsppi_q(d, q);

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
            ds_wsppi_residual(x, b, r);
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
            printf("    reached %e, on iter %i, delta %e\n", tol_norm_b, iter, delta_new);
            break;
        }

        // Output progress.
        printf("    iter %i, delta %e\n", iter, delta_new);
    }

    // don't forget to free the work arrays.
    delete [] d;
    delete [] r;
    delete [] q;
    delete [] s;
}

// ===
// Section for S^{mp} x
// ===

void
ds_smp_x(const double * const x, double * const m)
{
    const DSPoint * const l_mc = map_coords;
    const DSPixel * const l_pixels = pixels;

    const double var_s = p.corr_var_s;
    const double l_para_2 = p.corr_l_para_2;
    const double l_perp_2 = p.corr_l_perp_2;

    const int gt_n = gt.n;
    const double gt_dx = gt.dx;
    const double * const gt_table = gt.table;

    const int map_n = p.map_n;
    const int pix_n = p.num_pixels;
    int i, j;

#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < map_n; ++i) {
        m[i] = 0.0;
        double xi = l_mc[i].x;
        double yi = l_mc[i].y;
        double zi = l_mc[i].z;
        for (j = 0; j < pix_n; ++j) {
            double dx = xi - l_pixels[j].x;
            double dy = yi - l_pixels[j].y;
            double dz = zi - l_pixels[j].z;
            double x_perp_2 = (dx*dx + dy*dy) / l_perp_2;
            double x_para_2 = dz*dz / l_para_2;
            double g_perp = ds_linterp(x_perp_2, gt_n, gt_dx, gt_table);
            double g_para = ds_linterp(x_para_2, gt_n, gt_dx, gt_table);
            double smp_ij = var_s * g_perp * g_para;
            m[i] += smp_ij * x[j];
        }
    }
}

void
ds_compute_map(const double * const d, double * const m)
{
    const int pix_n = p.num_pixels;
    const int map_n = p.map_n;

    // Where we will store the pcg solution...
    double * const x = new double[pix_n];
    // Zero out for initial solution.
    for (int i = 0; i < pix_n; ++i) {
        x[i] = 0.0;
    }

    printf("Calling PCG with n_pix = %i\n", pix_n);
    Timer *t1 = new Timer();

    // x = (W S^{pp} + I)^{-1} W d
    pcg_wsppi_b(d, x);

    printf("[TIME] %g ms for x = (W S^{pp} + I)^{-1} W d.\n", t1->elapsed());
    delete t1;

    printf("Calling map multiply with map_n = %i, n = %.2e.\n", map_n, (double)map_n * pix_n);
    Timer *t2 = new Timer();

    // m = S^{mp} x
    ds_smp_x(x, m);

    printf("[TIME] %g ms for m = S^{mp} x.\n", t2->elapsed());
    delete t2;

    delete [] x;
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
