/*
2014, the dachshund authors.
*/

#ifndef __DS_PCG_H__
#define __DS_PCG_H__

/*
Linear interpolation, given x and table data.

We assume x >= 0. If x is larger than the last tabulated point, we return 0.

Parameters
x : the eval point
nx : the number of samples
dx: the sample point separation
f_table : array containing the function value at sample points.
*/
static inline double
ds_linterp(const double x, const int nx, const double dx,
    const double * const f_table)
{
    // in index space.
    double xx = x / dx;
    // decimate to get low edge.
    int ix = xx;
    // Check if we are outside table range.
    if (ix > nx - 2) {
        return 0.0;
    }
    // the linterp expression.
    return f_table[ix] * (ix + 1 - xx) + f_table[ix+1] * (xx - ix);
}

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

extern DSFuncTable gt;
extern DSPixel *pixels;
extern DSPoint *map_coords;


void
ds_lum_init(const double * const skewer_x, const double * const skewer_y,
    const double * const pixel_weights);

void
ds_lum_free();

void
pcg_wsppi(const int n, double * const x, const double * const b,
    const int max_iter, const double tol, const int verbose);


void
pcg_wsppi_spm(const int n, const int k, double * const x,
    const int max_iter, const double tol, const int verbose);

#endif
