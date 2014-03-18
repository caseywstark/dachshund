/*
2014, the dachshund authors.
*/

#include <cmath>
#include <cstdio>

#include "pcg.h"

void
pcg(const int n, double (*A)(int, int), double * const x,
    const double * const b, const int max_iter, const double tol,
    const int verbose)
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
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        double Ax_i = 0.0;
        for (j = 0; j < n; ++j) {
            Ax_i += (*A)(i, j) * x[j];
        }
        r[i] = b[i] - Ax_i;
    }

    // The preconditioning.
    // For now, we will just try a Jacobian M (A diag).
    // s = M^-1 r
    for (i = 0; i < n; ++i) {
        d[i] = r[i] / (*A)(i, i);
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
#if defined(_OPENMP)
        #pragma omp parallel for private(i, j)
#endif
        for (i = 0; i < n; ++i) {
            q[i] = 0.0;
            for (j = 0; j < n; ++j) {
                q[i] += (*A)(i, j) * d[j];
            }
        }

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
#if defined(_OPENMP)
            #pragma omp parallel for private(i, j)
#endif
            for (i = 0; i < n; ++i) {
                double Ax_i = 0.0;
                for (j = 0; j < n; ++j) {
                    Ax_i += (*A)(i, j) * x[j];
                }
                r[i] = b[i] - Ax_i;
            }
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
            s[i] = r[i] / (*A)(i, i);
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
pcg_covar(const int n, double (*A)(int, int), double * const x,
    double (*b)(int, int), const int bj, const int max_iter, const double tol,
    const int verbose)
{
    int i, j;

    // Work arrays
    double *d = new double[n];
    double *r = new double[n];
    double *s = new double[n];
    double *q = new double[n];

    // Setup the residual.
    // r = b - A x
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < n; ++i) {
        double Ax_i = 0.0;
        for (j = 0; j < n; ++j) {
            Ax_i += (*A)(i, j) * x[j];
        }
        r[i] = (*b)(i, bj) - Ax_i;
    }

    // The preconditioning.
    // For now, we will just try a Jacobian M (A diag).
    // s = M^-1 r
    for (i = 0; i < n; ++i) {
        d[i] = r[i] / (*A)(i, i);
    }

    // delta = r^T d
    double delta_new = 0.0;
    for (i = 0; i < n; ++i) {
        delta_new += r[i] * d[i];
    }

    if (verbose) {
        printf("  starting PCG with delta %e\n", delta_new);
    }

    for (int iter = 1; iter < max_iter; ++iter) {
        // q = A d
#if defined(_OPENMP)
        #pragma omp parallel for private(i, j)
#endif
        for (i = 0; i < n; ++i) {
            q[i] = 0.0;
            for (j = 0; j < n; ++j) {
                q[i] += (*A)(i, j) * d[j];
            }
        }

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
#if defined(_OPENMP)
            #pragma omp parallel for private(i, j)
#endif
            for (i = 0; i < n; ++i) {
                double Ax_i = 0.0;
                for (j = 0; j < n; ++j) {
                    Ax_i += (*A)(i, j) * x[j];
                }
                r[i] = (*b)(i, bj) - Ax_i;
            }
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
            s[i] = r[i] / (*A)(i, i);
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
        if (delta_new > tol) {
            if (verbose) {
                printf("  reached tolerance %e, on iter %i, delta %e\n", tol, iter, delta_new);
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
