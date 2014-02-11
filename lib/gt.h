/*
2014, the dachshund authors.

GT is short for gaussian table.
*/

#ifndef __DS_GT_H__
#define __DS_GT_H__

double
linterp(const double x, const int nx, const double dx,
    const double * const f_table);

class GT {
  public:
    GT(const int n, const double dx);
    ~GT();

    int n;
    double dx;
    double *table;

    double operator()(const double x);
};

#endif
