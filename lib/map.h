/*
2014, the dachshund authors.
*/

#ifndef __DS_MAP_H__
#define __DS_MAP_H__

#include <cmath>
#include <cstdio>

#include <string>
#include <vector>

#include "pixel.h"
#include "signal_covar.h"
#include "wf_vector.h"

class Map {
 public:
  const int n_pix;
  const Pixel* const pixels;
  const int n_map;
  const Point* const map_coords;
  const SignalCovarParams* const s_params;
  const WFVector* const wf_vector;

  Map(const int n_pix, const Pixel* const pixels,
      const int n_map, const Point* const map_coords,
      const SignalCovarParams* const s_params,
      const WFVector* const wf_vector)
      : n_pix(n_pix), pixels(pixels), n_map(n_map), map_coords(map_coords),
        s_params(s_params), wf_vector(wf_vector) {  }

  void set_wf_vector(double * const x) const;
  void set_map(double * const m, const double * const x) const;

  void solve(double * const m) const;
  void solve_covar_diag(double * const sigma_m) const;
};

void smp_product(const int num_pixel_points, const Point* const pixel_coords,
    const int num_map_points, const Point* const map_coords,
    const SignalCovarParams* const s_params,
    const double* const x, double* const m);

#endif
