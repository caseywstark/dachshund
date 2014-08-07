/*
2014, the dachshund authors.
*/

#ifndef __DS_UTIL_H__
#define __DS_UTIL_H__

#include <string>

#include "map.h"

void
read_pixel_data(const std::string path, const int num_pixels,
    Pixel * const pixels);

void
init_uniform_map_coords(const double lx, const double ly, const double lz,
    const int nx, const int ny, const int nz, Point * const coords);

void
check_finite(const int n, const double * const v);

#endif
