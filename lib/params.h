/*
2014, the dachshund authors.
*/

#ifndef __DS_PARAMS_H__
#define __DS_PARAMS_H__

#include <string>

#include "gt.h"

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

struct ds_params_s {
    // cells
    int nx, ny, nz;
    int num_cells;
    double lx, ly, lz;
    double dx, dy, dz;

    // skewers
    int num_skewers;
    int num_pixels;
    int num_pixel_elements;
    double dz_pix;

    // correlation
    double l_perp;
    double l_para;
    double sigma;

    // PCG
    int pcg_max_iter;
    double pcg_tol;

    // options
    int compute_covar;

    // data access
    DSPixel *pixels;
    DSPoint *map;

    // paths
    std::string skewer_x_path;
    std::string skewer_y_path;
    std::string pixel_data_path;
    std::string pixel_weights_path;
    std::string map_path;

    // misc
    GT *gt;
};

typedef struct ds_params_s DSParams;

extern DSParams p;

void
ds_params_init(const std::string config_path);

#endif
