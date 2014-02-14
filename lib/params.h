/*
2014, the dachshund authors.
*/

#ifndef __DS_PARAMS_H__
#define __DS_PARAMS_H__

#include "gt.h"

struct ds_params_s {
    // cells
    int nx, ny, nz;
    double lx, ly, lz;
    double dx, dy, dz;

    // skewers
    int num_skewers;
    int num_pixels;
    double dz_pix;

    // correlation
    double l_perp;
    double l_para;
    double sigma;

    // PCG
    int pcg_max_iter;
    double pcg_tol;

    // data access
    double *skewer_x;
    double *skewer_y;
    double *noise;

    // misc
    GT *gt;
};

typedef struct ds_params_s DSParams;

extern DSParams p;

#endif
