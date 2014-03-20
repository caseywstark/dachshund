/*
2014, the dachshund authors.
*/

#ifndef __DS_PCG_H__
#define __DS_PCG_H__

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

void
ds_lum_init(const double * const skewer_x, const double * const skewer_y,
    const double * const pixel_weights);

void
ds_lum_free();

void
pcg_wsppi(const int n, double * const x, const double * const b,
    const int max_iter, const double tol, const int verbose);

void
ds_smp_x(const double * const x, double * const m);

void
pcg_wsppi_spm(const int n, const int k, double * const x,
    const int max_iter, const double tol, const int verbose);

#endif
