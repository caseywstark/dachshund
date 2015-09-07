
#include <cmath>
#include <cstdio>

#include <iostream>
#include <Eigen/Dense>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "linalg.h"

#include "wf_vector.h"

double WFVectorSN::A(const int i, const int j, const bool exact) {
  const int delta_ij = (i == j);

  const double x_i = pixels[i].x;
  const double y_i = pixels[i].y;
  const double z_i = pixels[i].z;
  const double N_ii = pixels[i].N;

  const double dx = x_i - pixels[j].x;
  const double dy = y_i - pixels[j].y;
  const double dz = z_i - pixels[j].z;

  const double x_perp_2 = dx*dx + dy*dy;
  const double x_para_2 = dz*dz;
  double S_ij;
  if (exact) {
    S_ij = signal_covar_exact(x_perp_2, x_para_2, s_params);
  }
  else {
    S_ij = signal_covar(x_perp_2, x_para_2, s_params);
  }

  const double A_ij = S_ij + N_ii * delta_ij;
  return A_ij;
}

void WFVectorSN::preconditioner_solve(const double* const b, double* const x) {
  const double var_s = s_params->var_s;
  for (int i = 0; i < n; ++i) {
    // no coord lookup for diagonal. Apc_ii = var_s + N_i
    double Apc_ii = var_s + pixels[i].N;
    x[i] = b[i] / Apc_ii;
  }
}

void WFVectorSN::A_product(const double* const a, double* const b) {
  int i, j;
#if defined(_OPENMP)
  #pragma omp parallel for private(i, j)
#endif
  for (i = 0; i < n; ++i) {
    b[i] = 0.0;
    const double x_i = pixels[i].x;
    const double y_i = pixels[i].y;
    const double z_i = pixels[i].z;
    const double N_ii = pixels[i].N;
    for (j = 0; j < n; ++j) {
      const int delta_ij = (i == j);
      const double dx = x_i - pixels[j].x;
      const double dy = y_i - pixels[j].y;
      const double dz = z_i - pixels[j].z;
      const double x_perp_2 = (dx*dx + dy*dy);
      const double x_para_2 = dz*dz;
      const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
      const double A_ij = S_ij + N_ii * delta_ij;
      b[i] += A_ij * a[j];
    }
  }
}

void WFVectorSN::A_residual(const double* const x, double* const r) {
  int i, j;
#if defined(_OPENMP)
  #pragma omp parallel for private(i, j)
#endif
  for (i = 0; i < n; ++i) {
    double sum_Ax_i = 0.0;
    const double x_i = pixels[i].x;
    const double y_i = pixels[i].y;
    const double z_i = pixels[i].z;
    const double N_ii = pixels[i].N;
    const double d_i = pixels[i].d;

    for (j = 0; j < n; ++j) {
      const int delta_ij = (i == j);
      const double dx = x_i - pixels[j].x;
      const double dy = y_i - pixels[j].y;
      const double dz = z_i - pixels[j].z;
      const double x_perp_2 = (dx*dx + dy*dy);
      const double x_para_2 = dz*dz;
      const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
      const double A_ij = S_ij + N_ii * delta_ij;
      sum_Ax_i += A_ij * x[j];
    }
    // r_i = b_i - sum_j A_ij x_j
    r[i] = d_i - sum_Ax_i;
  }
}

double WFVectorSN::norm_b() {
  double sum_b2 = 0.0;
  for (int i = 0; i < n; ++i) {
    const double bi = pixels[i].d;
    sum_b2 += bi * bi;
  }
  return sqrt(sum_b2);
}

PCGResult WFVectorSN::solve_pcg(const PCGParams& pcg_params,
    double* const x) {
  const double norm_b_val = norm_b();
  PCGResult res = pcg_solve(pcg_params, norm_b_val, this_A_product,
      this_A_residual, this_pc_solve, x);

  return res;
}

void WFVectorSN::solve_cf(double* const x, const bool verbose, const bool exact) {
  using namespace Eigen;
  if (verbose) { puts("Filling A matrix."); }

  // init A = (S + N)
  MatrixXd sn(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      sn(i, j) = A(i, j, exact);
    }
  }

  if (verbose) { puts("Filling b vector."); }

  // init b = d
  VectorXd b(n);
  for (int i = 0; i < n; ++i) {
    b(i) = pixels[i].d;
  }

  if (verbose) { puts("Calling eigen solve."); }

  // solve
  VectorXd vec_x(n);
  vec_x = sn.ldlt().solve(b);

  if (verbose) { puts("Copying back."); }

  // write back to our format.
  for (int i = 0; i < n; ++i) {
    x[i] = vec_x(i);
  }
}

double WFVectorW::A(const int i, const int j, const bool exact) {
  const int delta_ij = (i == j);

  const double x_i = pixels[i].x;
  const double y_i = pixels[i].y;
  const double z_i = pixels[i].z;
  const double w_i = pixels[i].w;

  const double dx = x_i - pixels[j].x;
  const double dy = y_i - pixels[j].y;
  const double dz = z_i - pixels[j].z;
  const double w_j = pixels[j].w;

  const double x_perp_2 = (dx*dx + dy*dy);
  const double x_para_2 = dz*dz;

  double S_ij;
  if (exact) {
    S_ij = signal_covar_exact(x_perp_2, x_para_2, s_params);
  }
  else {
    S_ij = signal_covar(x_perp_2, x_para_2, s_params);
  }

  const double A_ij = w_i * S_ij * w_j + delta_ij;
  return A_ij;
}

void WFVectorW::preconditioner_solve(const double* const b, double* const x) {
  const double var_s = s_params->var_s;

  for (int i = 0; i < n; ++i) {
    // no coord lookup for diagonal. Apc_ii = var_s + N_i
    double wi = pixels[i].w;
    double Apc_ii = wi * wi * var_s + 1.0;
    x[i] = b[i] / Apc_ii;
  }
}

void WFVectorW::A_product(const double* const a, double* const b) {
  int i, j;
#if defined(_OPENMP)
  #pragma omp parallel for private(i, j)
#endif
  for (i = 0; i < n; ++i) {
    b[i] = 0.0;
    const double x_i = pixels[i].x;
    const double y_i = pixels[i].y;
    const double z_i = pixels[i].z;
    const double w_i = pixels[i].w;
    for (j = 0; j < n; ++j) {
      const int delta_ij = (i == j);
      const double w_j = pixels[j].w;
      const double dx = x_i - pixels[j].x;
      const double dy = y_i - pixels[j].y;
      const double dz = z_i - pixels[j].z;
      const double x_perp_2 = (dx*dx + dy*dy);
      const double x_para_2 = dz*dz;
      const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
      const double A_ij = w_i * S_ij * w_j + delta_ij;
      b[i] += A_ij * a[j];
    }
  }
}

void WFVectorW::A_residual(const double* const x, double* const r) {
  int i, j;
#if defined(_OPENMP)
  #pragma omp parallel for private(i, j)
#endif
  for (i = 0; i < n; ++i) {
    double sum_Ax_i = 0.0;
    const double x_i = pixels[i].x;
    const double y_i = pixels[i].y;
    const double z_i = pixels[i].z;
    const double w_i = pixels[i].w;

    for (j = 0; j < n; ++j) {
      const int delta_ij = (i == j);
      const double w_j = pixels[j].w;
      const double dx = x_i - pixels[j].x;
      const double dy = y_i - pixels[j].y;
      const double dz = z_i - pixels[j].z;
      const double x_perp_2 = (dx*dx + dy*dy);
      const double x_para_2 = dz*dz;
      const double S_ij = signal_covar(x_perp_2, x_para_2, s_params);
      const double A_ij = w_i * S_ij * w_j + delta_ij;
      sum_Ax_i += A_ij * x[j];
    }

    // r_i = b_i - sum_j A_ij x_j
    const double wd_i = pixels[i].wd;
    r[i] = wd_i - sum_Ax_i;
  }
}

double WFVectorW::norm_b() {
  double sum_b2 = 0.0;
  for (int i = 0; i < n; ++i) {
    const double bi = pixels[i].wd;
    sum_b2 += bi * bi;
  }
  return sqrt(sum_b2);
}


PCGResult WFVectorW::solve_pcg(const PCGParams& pcg_params, double* const x) {
  const double norm_b_val = norm_b();
  PCGResult res = pcg_solve(pcg_params, norm_b_val, this_A_product,
      this_A_residual, this_pc_solve, x);

  // Don't forget the last multiply!
  for (int i = 0; i < n; ++i) {
    x[i] *= pixels[i].w;
  }

  return res;
}

void WFVectorW::solve_cf(double* const x, const bool verbose,
    const bool exact) {
  using namespace Eigen;

  if (verbose) {
    puts("Filling A matrix.");
  }

  // init A = (S + N)
  MatrixXd sn(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      sn(i, j) = A(i, j, exact);
    }
  }

  if (verbose) {
    puts("Filling b vector.");
  }

  // init b = d
  VectorXd b(n);
  for (int i = 0; i < n; ++i) {
    b(i) = pixels[i].wd;
  }

  if (verbose) {
    puts("Calling eigen solve.");
  }

  // solve
  VectorXd vec_x(n);
  vec_x = sn.ldlt().solve(b);

  if (verbose) {
    puts("Copying back.");
  }

  // write back to our format
  // Don't forget the last multiply!
  for (int i = 0; i < n; ++i) {
    x[i] = pixels[i].w * vec_x(i);
  }
}
