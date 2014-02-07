/*
2014, the dachshund authors.

GPT is short for gaussian product table.
*/

#ifndef __DS_GPT_H__
#define __DS_GPT_H__

double
bilinterp(const double x, const double y, const int nx, const int ny,
    const double dx, const double dy, const double * const f_table);

class GPT {
  public:
    GPT(const int n, const double dx);
    ~GPT();

    int n;
    double dx;
    double *table;

    double operator()(const double x, const double y);
};

#endif
