/*
2014, the dachshund authors.
*/

#include <cmath>
#include <cstdio>

#include <iostream>
#include <Eigen/Dense>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "timer.h"
#include "linalg.h"

#include "map.h"

void smp_product(const int num_pixel_points, const Point* const pixel_coords,
    const int num_map_points, const Point* const map_coords,
    const SignalCovarParams* const s_params,
    const double* const x, double* const m) {
  int i, j;
#if defined(_OPENMP)
  #pragma omp parallel for private(i, j)
#endif
  for (i = 0; i < num_map_points; ++i) {
    m[i] = 0.0;
    const double x_i = map_coords[i].x;
    const double y_i = map_coords[i].y;
    const double z_i = map_coords[i].z;
    for (j = 0; j < num_pixel_points; ++j) {
      const double dx = x_i - pixel_coords[j].x;
      const double dy = y_i - pixel_coords[j].y;
      const double dz = z_i - pixel_coords[j].z;
      const double x_perp_2 = (dx*dx + dy*dy);
      const double x_para_2 = dz*dz;
      const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
      m[i] += S_ij * x[j];
    }
  }
}
