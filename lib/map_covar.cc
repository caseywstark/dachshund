
#include "map_covar.h"

#include "wf_vector.h"

// Section for \sigma_m = diag(M) = diag(S^{mp} W S^{pm})
// This is the diagonal of the map covariance under the assumption N >> S.
// Since W is diagonal,
// \sigma_m_i = M_ii = S^mp_ij W_jk S^pm_kl = S^mp_ij W_jj S^pm_ji

void sigma_m_noise_dom(const int n_pix, const NPixel* pixels,
    const int n_map, const Point* map_coords,
    const SignalCovarParams* s_params, double* sigma_m) {
  int i, j;
#if defined(_OPENMP)
  #pragma omp parallel for private(i, j)
#endif
  for (i = 0; i < n_map; ++i) {
    // Reset value at this map point.
    sigma_m[i] = 0.0;
    // Grab map coords.
    const double xi = map_coords[i].x;
    const double yi = map_coords[i].y;
    const double zi = map_coords[i].z;

    for (j = 0; j < n_pix; ++j) {
      const double dx = xi - pixels[j].x;
      const double dy = yi - pixels[j].y;
      const double dz = zi - pixels[j].z;
      const double W_jj = 1.0 / pixels[j].N;
      const double x_perp_2 = dx*dx + dy*dy;
      const double x_para_2 = dz*dz;
      const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
      // since S_ij^mp = S_ji^pm
      sigma_m[i] += S_ij * W_jj * S_ij;
    }
  }
}

PCGResult MapCovarDiag::solve_column_pcg(const int i,
    const PCGParams& pcg_params, double* const x) {
  // setup 'pixels'
  const double xi = map_coords[i].x;
  const double yi = map_coords[i].y;
  const double zi = map_coords[i].z;

  NPixel *col_pixels = new NPixel[n_pix];
  for (int k = 0; k < n_pix; ++k) {
    const double dx = xi - pixels[k].x;
    const double dy = yi - pixels[k].y;
    const double dz = zi - pixels[k].z;
    const double x_perp_2 = dx*dx + dy*dy;
    const double x_para_2 = dz*dz;
    const double S_ik = signal_covar(x_perp_2, x_para_2, s_params);

    col_pixels[k].x = pixels[k].x;
    col_pixels[k].y = pixels[k].y;
    col_pixels[k].z = pixels[k].z;
    col_pixels[k].N = pixels[k].N;
    col_pixels[k].d = S_ik;
  }

  WFVectorSN *wf = new WFVectorSN(n_pix, col_pixels, s_params);

  PCGResult res = wf->solve_pcg(pcg_params, x);

  delete wf;
  delete [] col_pixels;

  return res;
}

std::vector<PCGResult> MapCovarDiag::solve_pcg(const PCGParams& pcg_params,
    double* const sigma_m) {
  std::vector<PCGResult> results;
  double *x = new double[n_pix];

  // Iterate over map points.
  for (int i = 0; i < n_map; ++i) {
    // Print progress because this part is sloooow.
    printf("cell %i -- %f.\n", i, (double)i/n_map);

    // Reset.
    sigma_m[i] = 0.0;
    for (int j = 0; j < n_pix; ++j) { x[j] = 0.0; }

    // Solve WFVector.
    PCGResult res = solve_column_pcg(i, pcg_params, x);
    results.push_back(res);

    // Annoying edge case...
    // Sometimes [Sum_k Spm_ki^2]^{1/2} = 0, meaning that the pixel has no
    // overlap with the map points.
    // In this case, M_ii = 0 anyway, so don't continue.
    if (res.success) {
      // Map multiply.
      const double xi = map_coords[i].x;
      const double yi = map_coords[i].y;
      const double zi = map_coords[i].z;
      for (int j = 0; j < n_pix; ++j) {
        const double dx = xi - pixels[j].x;
        const double dy = yi - pixels[j].y;
        const double dz = zi - pixels[j].z;
        const double x_perp_2 = dx*dx + dy*dy;
        const double x_para_2 = dz*dz;
        const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
        sigma_m[i] += S_ij * x[j];
      }
    }
  }

  return results;
}
