/*

*/

#ifndef __DS_LINALG_H__
#define __DS_LINALG_H__

double
ds_vector_dot(const int n, const double * const a, const double * const b);

// Could just call vector dot for this. Not sure if it's more optimized with
// only one pointer (aliasing stuff).
double
ds_vector_norm(const int n, const double * const x);

// c_i = c_i + a + b_i
void
ds_vector_add_cv(const int n, const double a, const double * const b,
    double * const c);

// c_i = c_i - a + b_i
void
ds_vector_sub_cv(const int n, const double a, const double * const b,
    double * const c);

#endif
