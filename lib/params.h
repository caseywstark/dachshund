/*
2014, the dachshund authors.
*/

#ifndef __DS_PARAMS_H__
#define __DS_PARAMS_H__

#include <string>

struct ds_params_s {
    // domain
    double lx, ly, lz;

    // skewers
    int num_skewers;
    int num_pixels;
    int pix_n;
    double pix_dz;

    // map
    int map_nx, map_ny, map_nz;
    int map_n;
    double map_dx, map_dy, map_dz;

    // correlation
    double corr_var_s;
    double corr_l_perp_2;
    double corr_l_para_2;

    // PCG
    int pcg_max_iter;
    double pcg_tol;

    // options
    int option_compute_covar;
    int option_smooth_map;

    // paths
    std::string skewer_x_path;
    std::string skewer_y_path;
    std::string pixel_data_path;
    std::string pixel_weights_path;
    std::string map_path;
};

typedef struct ds_params_s DSParams;

extern DSParams p;

void
ds_params_init(const std::string config_path);

#endif
