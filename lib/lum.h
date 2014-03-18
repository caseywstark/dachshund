/*
2014, the dachshund authors.
*/

#ifndef __DS_LUM_H__
#define __DS_LUM_H__

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
ds_lum_gt_init(const int n, const double dx);

void
ds_lum_gt_free();

void
ds_lum_pixels_init(const int num_skewers,
    const double * const skewer_x, const double * const skewer_y,
    const int num_pixels, const double pix_dz,
    const double * const pixel_weights);

void
ds_lum_pixels_free();

void
ds_lum_map_coords_init(const int map_nx, const int map_ny, const int map_nz,
    const int map_dx, const int map_dy, const int map_dz);

void
ds_lum_map_coords_free();

//
// Matrix lookup functions
//

double
smp_lookup(const int i, const int j);

double
wsppi_lookup(const int i, const int j);

double
wspm_lookup(const int i, const int j);

#endif
