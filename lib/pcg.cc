/*
2014, the dachshund authors.
*/

#include <math.h>
#include <stdint.h>
#include <cstdio>
#include <cstdlib>

#include <omp.h>

#include "lum.h"
#include "util.h"
#include "pcg.h"

/*
pcg_invert solves for the i-th component of A^-1.
*/
void
pcg_invert(const int k, LUM *A, double * const x,
    const int max_iter, const double tol,
    double * const d, double * const r,
    double * const s, double * const q)
{
    // Notation:
    //   k is the row we are inverting.
    //   A is the n x m matrix.
    //   x = A^{-1} is the solution.
    //   r is the residual.
    //   ...

    const int n = A->n;
    const int m = A->m;

    // Setup the residual.
    // r = b - A x

    int i, j;
    #pragma omp parallel for default(none) private(tid, i, j)
    for (i = 0; i < n; ++i) {
        // Handle b[i]
        if (i == k) {
            r[i] = 1.0;
        }
        else {
            r[i] = 0.0;
        }

        //r[i] -= A(i, j) * x[j];
        for (j = 0; j < m; ++j) {
            r[i] -= (*A)(i, j) * x[j];
        }
    }

    // The preconditioning.
    // For now, we will just try a Jacobian M (A diag).
    // s = M^-1 r
    for (int j = 0; j < m; ++j) {
        d[j] = r[j] / (*A)(j, j);
    }

    // delta = r^T d
    double delta_new = 0.0;
    for (int j = 0; j < m; ++j) {
        delta_new += r[j] * d[j];
    }

    double delta_0 = delta_new;
    printf("  delta_0 = %g\n", delta_0);

    // Handle singular case.
    if (delta_0 == 0.0) {
        return;
    }

    // the iteration count.
    int iter = 0;
    // the CG loop...
    while (iter < max_iter && delta_new > tol * delta_0) {
        // q = A d
        for (int i = 0; i < n; ++i) {
            q[i] = 0.0;
            for (int j = 0; j < m; ++j) {
                q[i] += (*A)(i, j) * d[j];
            }
        }

        // alpha = delta_new / (d^T q)
        double denom = 0.0;
        for (int i = 0; i < n; ++i) {
            denom += d[i] * q[i];
        }
        double alpha = delta_new / denom;

        // x = x + alpha d
        for (int j = 0; j < m; ++j) {
            x[j] += alpha * d[j];
        }

        // Update residual.
        if (iter != 0 && iter % 100 == 0) {
            // Full update.
            // r = b - A x
            for (int i = 0; i < n; ++i) {
                // Handle b[i]
                if (i == k) {
                    r[i] = 1.0;
                }
                else {
                    r[i] = 0.0;
                }

                //r[i] -= A(i, j) * x[j];
                for (int j = 0; j < m; ++j) {
                    r[i] -= (*A)(i, j) * x[j];
                }
            }
        }
        else {
            // Approximate update.
            // r = r - alpha q
            for (int i = 0; i < n; ++i) {
                r[i] -= alpha * q[i];
            }
        }

        // reapply preconditioner.
        for (int j = 0; j < m; ++j) {
            s[j] = r[j] / (*A)(j, j);
        }

        // save current delta.
        double delta_old = delta_new;

        // Update delta.
        // delta_new = r^T s
        delta_new = 0.0;
        for (int j = 0; j < m; ++j) {
            delta_new += r[j] * s[j];
        }

        // Update d.
        double beta = delta_new / delta_old;
        for (int j = 0; j < m; ++j) {
            d[j] = s[j] + beta * d[j];
        }

        //printf("    iter %i   delta %g\n", iter, delta_new);

        // Finally, update the count.
        ++iter;
    }

    // Handle case of exceeding max_iter.
    if (iter == max_iter && delta_new > tol * delta_0) {
        fprintf(stderr, "PCG convergence: did not reach solution in allowed iterations.\n");
        fprintf(stderr, "max_iter %i, tol %e, delta_0 %e, delta %e\n",
            max_iter, tol, delta_0, delta_new);
        exit(1);
    }

    // DEBUG
    printf("  PCG in %i iterations to delta / delta_0 %e.\n",
        iter, delta_new / delta_0);
}
