/*
2014, the dachshund authors.
*/

#include <cmath>
#include <cstdio>

#include "params.h"
#include "pcg.h"

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
ds_pixels_init(const int num_skewers,
    const double * const skewer_x, const double * const skewer_y,
    const int num_pixels, const double pix_dz,
    const double * const pixel_weights)
{
    int pix_n = num_skewers * num_pixels;
    pixels = new DSPixel[pix_n];

    for (int i = 0; i < pix_n; ++i) {
        int isk = i / num_pixels;
        int iz = i % num_pixels;
        pixels[i].x = skewer_x[isk];
        pixels[i].y = skewer_y[isk];
        pixels[i].z = pix_dz * (iz + 0.5);
        pixels[i].w = pixel_weights[i];
    }
}

void
ds_pixels_free()
{
    delete pixels;
}

void
ds_map_coords_init(const int map_nx, const int map_ny, const int map_nz,
    const int map_dx, const int map_dy, const int map_dz)
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
ds_lum_init(const double * const skewer_x, const double * const skewer_y,
    const double * const pixel_weights)
{
    ds_pixels_init(p.num_skewers, skewer_x, skewer_y, p.num_pixels,
        p.pix_dz, pixel_weights);
    ds_map_coords_init(p.map_nx, p.map_ny, p.map_nz, p.map_dx, p.map_dy,
        p.map_dz);

    int gt_n = 50;
    double gt_dx = 2.0 / gt_n;
    ds_gt_init(gt_n, gt_dx);
}

void
ds_lum_free()
{
    ds_map_coords_free();
    ds_pixels_free();
    ds_gt_free();
}

void
ds_wsppi_residual(const int n, const double * const x,
    const double * const b, double * const r)
{
    int i, j;
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        double Ax_i = 0.0;

        double xi = pixels[i].x;
        double yi = pixels[i].y;
        double zi = pixels[i].z;
        double one_minus_wi = 1.0 - pixels[i].w;

        for (j = 0; j < n; ++j) {
            // Separation vector components.
            double dx = xi - pixels[j].x;
            double dy = yi - pixels[j].y;
            double dz = zi - pixels[j].z;

            // The gaussian terms.
            double x_perp_2 = (dx*dx + dy*dy) / p.corr_l_perp_2;
            double x_para_2 = dz*dz / p.corr_l_para_2;

            // S = sigma^2 exp(...) exp(...)
            double g_perp = ds_linterp(x_perp_2, gt.n, gt.dx, gt.table);
            double g_para = ds_linterp(x_para_2, gt.n, gt.dx, gt.table);
            double wsppi_ij = p.corr_var_s * g_perp * g_para;

            // nice trick to avoid branching.
            int delta_ij = (i == j);
            wsppi_ij = wsppi_ij * (1.0 - one_minus_wi * delta_ij) + delta_ij;

            // Update sum.
            Ax_i += wsppi_ij * x[j];
        }
        r[i] = b[i] - Ax_i;
    }
}

void
ds_wsppi_q(const int n, const double * const d, double * const q)
{
    int i, j;
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        q[i] = 0.0;
        double xi = pixels[i].x;
        double yi = pixels[i].y;
        double zi = pixels[i].z;
        double one_minus_wi = 1.0 - pixels[i].w;
        for (j = 0; j < n; ++j) {
            double dx = xi - pixels[j].x;
            double dy = yi - pixels[j].y;
            double dz = zi - pixels[j].z;
            double x_perp_2 = (dx*dx + dy*dy) / p.corr_l_perp_2;
            double x_para_2 = dz*dz / p.corr_l_para_2;
            double g_perp = ds_linterp(x_perp_2, gt.n, gt.dx, gt.table);
            double g_para = ds_linterp(x_para_2, gt.n, gt.dx, gt.table);
            double wsppi_ij = p.corr_var_s * g_perp * g_para;
            int delta_ij = (i == j);
            wsppi_ij = wsppi_ij * (1.0 - one_minus_wi * delta_ij) + delta_ij;
            q[i] += wsppi_ij * d[j];
        }
    }
}

void
pcg_wsppi(const int n, double * const x, const double * const b,
    const int max_iter, const double tol, const int verbose)
{
    int i, j;

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
    ds_wsppi_residual(n, x, b, r);

    // The preconditioning.
    // For now, we will just try a Jacobian M (A diag).
    // s = M^-1 r
    for (i = 0; i < n; ++i) {
        // no coord lookup for diagonal piece...
        // just w_i var + 1
        double wsppi_ii = pixels[i].w * p.corr_var_s + 1.0;
        d[i] = r[i] / wsppi_ii;
    }

    // delta = r^T d
    double delta_new = 0.0;
    for (i = 0; i < n; ++i) {
        delta_new += r[i] * d[i];
    }

    if (verbose) {
        printf("  starting PCG with delta %e, tol * |b| %e\n", delta_new, tol_norm_b);
    }

    for (int iter = 1; iter < max_iter; ++iter) {
        // q = A d
        ds_wsppi_q(n, d, q);

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
        if (iter % 10 == 0) {
            // Full update.
            // r = b - A x
            ds_wsppi_residual(n, x, b, r);
        }
        else {
            // Approximate update.
            // r = r - alpha q
            for (j = 0; j < n; ++j) {
                r[j] -= alpha * q[j];
            }
        }

        // reapply preconditioner.
        for (i = 0; i < n; ++i) {
            double wsppi_ii = pixels[i].w * p.corr_var_s + 1.0;
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
            if (verbose) {
                printf("  reached %e, on iter %i, delta %e\n", tol_norm_b, iter, delta_new);
            }
            break;
        }

        // Output progress.
        if (verbose) {
            printf("  iter %i, delta %e\n", iter, delta_new);
        }
    }

    // don't forget to free the work arrays.
    delete [] d;
    delete [] r;
    delete [] q;
    delete [] s;
}

void
ds_wsppi_spm_residual(const int n, const int k, const double * const x,
    double * const r)
{
    int i, j;
    double xk = map_coords[k].x;
    double yk = map_coords[k].y;
    double zk = map_coords[k].z;

#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        double Ax_i = 0.0;
        double xi = pixels[i].x;
        double yi = pixels[i].y;
        double zi = pixels[i].z;
        double one_minus_wi = 1.0 - pixels[i].w;
        for (j = 0; j < n; ++j) {
            double dx = xi - pixels[j].x;
            double dy = yi - pixels[j].y;
            double dz = zi - pixels[j].z;
            double x_perp_2 = (dx*dx + dy*dy) / p.corr_l_perp_2;
            double x_para_2 = dz*dz / p.corr_l_para_2;
            double g_perp = ds_linterp(x_perp_2, gt.n, gt.dx, gt.table);
            double g_para = ds_linterp(x_para_2, gt.n, gt.dx, gt.table);
            double wsppi_ij = p.corr_var_s * g_perp * g_para;
            int delta_ij = (i == j);
            wsppi_ij = wsppi_ij * (1.0 - one_minus_wi * delta_ij) + delta_ij;
            Ax_i += wsppi_ij * x[j];
        }
        double dx = xi - xk;
        double dy = yi - yk;
        double dz = zi - zk;
        double x_perp_2 = (dx*dx + dy*dy) / p.corr_l_perp_2;
        double x_para_2 = dz*dz / p.corr_l_para_2;
        double g_perp = ds_linterp(x_perp_2, gt.n, gt.dx, gt.table);
        double g_para = ds_linterp(x_para_2, gt.n, gt.dx, gt.table);
        double spm_ik = p.corr_var_s * g_perp * g_para;

        r[i] = spm_ik - Ax_i;
    }
}

void
ds_wsppi_spm_q(const int n, const double * const d, double * const q)
{
    int i, j;
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        q[i] = 0.0;
        double xi = pixels[i].x;
        double yi = pixels[i].y;
        double zi = pixels[i].z;
        double one_minus_wi = 1.0 - pixels[i].w;
        for (j = 0; j < n; ++j) {
            double dx = xi - pixels[j].x;
            double dy = yi - pixels[j].y;
            double dz = zi - pixels[j].z;
            double x_perp_2 = (dx*dx + dy*dy) / p.corr_l_perp_2;
            double x_para_2 = dz*dz / p.corr_l_para_2;
            double g_perp = ds_linterp(x_perp_2, gt.n, gt.dx, gt.table);
            double g_para = ds_linterp(x_para_2, gt.n, gt.dx, gt.table);
            double wsppi_ij = p.corr_var_s * g_perp * g_para;
            int delta_ij = (i == j);
            wsppi_ij = wsppi_ij * (1.0 - one_minus_wi * delta_ij) + delta_ij;
            q[i] += wsppi_ij * d[j];
        }
    }
}

void
pcg_wsppi_spm(const int n, const int k, double * const x,
    const int max_iter, const double tol, const int verbose)
{
    int i, j;

    // Work arrays
    double *d = new double[n];
    double *r = new double[n];
    double *s = new double[n];
    double *q = new double[n];

    // Compute norm b for stop condition.
    double norm_b = 0.0;
    double xk = map_coords[k].x;
    double yk = map_coords[k].y;
    double zk = map_coords[k].z;

#if defined(_OPENMP)
    #pragma omp parallel for private(i)
#endif
    for (i = 0; i < n; ++i) {
        double xi = pixels[i].x;
        double yi = pixels[i].y;
        double zi = pixels[i].z;
        double dx = xi - xk;
        double dy = yi - yk;
        double dz = zi - zk;
        double x_perp_2 = (dx*dx + dy*dy) / p.corr_l_perp_2;
        double x_para_2 = dz*dz / p.corr_l_para_2;
        double g_perp = ds_linterp(x_perp_2, gt.n, gt.dx, gt.table);
        double g_para = ds_linterp(x_para_2, gt.n, gt.dx, gt.table);
        double spm_ik = p.corr_var_s * g_perp * g_para;
        norm_b += spm_ik * spm_ik;
    }

    norm_b = sqrt(norm_b);
    double tol_norm_b = tol * norm_b;

    // edge case...
    if (norm_b == 0.0) {
        return;
    }

    // Setup the residual.
    // r = b - A x
    ds_wsppi_spm_residual(n, k, x, r);

    // The preconditioning.
    // For now, we will just try a Jacobian M (A diag).
    // s = M^-1 r
    for (i = 0; i < n; ++i) {
        double wsppi_i = pixels[i].w * p.corr_var_s + 1.0;
        d[i] = r[i] / wsppi_i;
    }

    // delta = r^T d
    double delta_new = 0.0;
    for (i = 0; i < n; ++i) {
        delta_new += r[i] * d[i];
    }

    if (verbose) {
        printf("  starting PCG with delta %e, tol %e\n", delta_new, tol_norm_b);
    }

    for (int iter = 1; iter < max_iter; ++iter) {
        // q = A d
        ds_wsppi_spm_q(n, d, q);

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
        if (iter % 10 == 0) {
            // Full update.
            // r = b - A x
            ds_wsppi_spm_residual(n, k, x, r);
        }
        else {
            // Approximate update.
            // r = r - alpha q
            for (j = 0; j < n; ++j) {
                r[j] -= alpha * q[j];
            }
        }

        // reapply preconditioner.
        for (i = 0; i < n; ++i) {
            double wsppi_ii = pixels[i].w * p.corr_var_s + 1.0;
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
            if (verbose) {
                printf("  reached tolerance %e, on iter %i, delta %e\n", tol_norm_b, iter, delta_new);
            }
            break;
        }

        // Output progress.
        if (verbose) {
            printf("  iter %i, delta %e\n", iter, delta_new);
        }
    }

    // don't forget to free the work arrays.
    delete [] d;
    delete [] r;
    delete [] q;
    delete [] s;
}
