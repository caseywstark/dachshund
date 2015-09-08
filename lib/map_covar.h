
#ifndef __DS_MAP_COVAR_H__
#define __DS_MAP_COVAR_H__

#include <vector>

#include "pcg.h"
#include "pixel.h"
#include "signal_covar.h"

void sigma_m_noise_dom(const int n_pix, const NPixel* pixels,
    const int n_map, const Point* map_coords,
    const SignalCovarParams* s_params, double* sigma_m);

// The map covariance is
// M = S^mp (S^pp + N) S^pm .
// The map covariance diagonal is
// \sigma_m = diag(M) = S^mp (S^pp + N)^{-1} S^pm
// or
// \sigma_m_i = S^mp_ij (S^pp_jk + N_jk) S^pm_ki .

class MapCovarDiag {
 public:
  const int n_pix;
  const NPixel* const pixels;
  const int n_map;
  const Point* const map_coords;
  const SignalCovarParams* const s_params;

  MapCovarDiag(const int n_pix, const NPixel* const pixels,
      const int n_map, const Point* const map_coords,
      const SignalCovarParams* const s_params)
      : n_pix(n_pix), pixels(pixels), n_map(n_map), map_coords(map_coords),
        s_params(s_params) {  }

  PCGResult solve_column_pcg(const int i, const PCGParams& pcg_params,
      double* const x);
  std::vector<PCGResult> solve_pcg(const PCGParams& pcg_params,
      double* const sigma_m);
};

#endif
