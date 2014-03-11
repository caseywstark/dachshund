/*
2014, the dachshund authors.
*/

#ifndef __DS_PCG_H__
#define __DS_PCG_H__

void
pcg(const int n, double (*A)(int, int), double * const x,
    const double * const b, const int max_iter, const double tol);

void
pcg_covar(const int n, double (*A)(int, int), double * const x,
    double (*b)(int, int), const int bj, const int max_iter, const double tol);

#endif
