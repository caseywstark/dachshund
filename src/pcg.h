/*
2014, the dachshund authors.
*/

#ifndef __DS_PCG_H__
#define __DS_PCG_H__

void
pcg_invert(const int k, LUM *A, double * const x,
    const int max_iter, const double tol,
    double * const d, double * const r,
    double * const s, double * const q);

#endif
