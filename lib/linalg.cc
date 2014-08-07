/*

*/

#include "linalg.h"

#include <math.h>

double
ds_vector_dot(const int n, const double * const a, const double * const b)
{
    double dot = 0.0;
    for (int i = 0; i < n; ++i) {
        dot += a[i] * b[i];
    }
    return dot;
}

// Could just call vector dot for this. Not sure if it's more optimized with
// only one pointer (aliasing stuff).
double
ds_vector_norm(const int n, const double * const x)
{
    double sum_x2 = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_x2 += x[i]*x[i];
    }
    return sqrt(sum_x2);
}

// c_i = c_i + a + b_i
void
ds_vector_add_cv(const int n, const double a, const double * const b,
    double * const c)
{
    for (int i = 0; i < n; ++i) {
        c[i] += a * b[i];
    }
}

// c_i = c_i - a + b_i
void
ds_vector_sub_cv(const int n, const double a, const double * const b,
    double * const c)
{
    for (int i = 0; i < n; ++i) {
        c[i] -= a * b[i];
    }
}


