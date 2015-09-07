
#ifndef __DS_PCG_H__
#define __DS_PCG_H__

#include <cmath>
#include <functional>
#include <vector>

#include "linalg.h"

// PCG data structs

struct PCGParams {
  int n;
  int max_iter;
  int step_r;
  double tol;
  bool verbose;
};

struct PCGResult {
  bool success;
  int num_iters;
  std::vector<double> residual_norms;
};

static void pcg_step(const int n,
    std::function<void(const double*, double*)> A_product,
    std::function<void(const double*, double*)> A_residual,
    std::function<void(const double*, double*)> preconditioner_solve,
    double* const d, double* const q, double* const r,
    double* const s, double* const x, const bool full_residual,
    double* const delta) {
  // q = A d
  A_product(d, q);
  // alpha = delta / (d^T q)
  double alpha = *delta / ds_vector_dot(n, d, q);
  // Solution update. x = x + alpha d
  ds_vector_add_cv(n, alpha, d, x);

  // Update residual.
  if (full_residual) {
    // Full update, r = b - A x
    A_residual(x, r);
  }
  else {
    // Approximate update, r = r - alpha q
    ds_vector_sub_cv(n, alpha, q, r);
  }

  // reapply preconditioner, s = M^-1 r
  preconditioner_solve(r, s);

  // save current delta.
  double delta_old = *delta;

  // Update delta, delta_new = r^T s
  *delta = ds_vector_dot(n, r, s);

  // Update d, d = s + beta d
  double beta = *delta / delta_old;
  for (int i = 0; i < n; ++i) {
    d[i] = s[i] + beta * d[i];
  }
}

static PCGResult pcg_solve(const PCGParams& params, const double norm_b,
    std::function<void(const double*, double*)> A_product,
    std::function<void(const double*, double*)> A_residual,
    std::function<void(const double*, double*)> preconditioner_solve,
    double* const x) {
  const int n = params.n;
  const int max_iter = params.max_iter;
  const double tol = params.tol;
  const int step_r = params.step_r;
  const bool verbose = params.verbose;

  // return val.
  PCGResult res;

  // Catch special case.
  if (norm_b <= 0.0) {
    res.success = false;
    puts("[WARNING] Requested PCG solve where |b| <= 0.0. Returning early.\n");
    return res;
  }

  // Work arrays
  double *d = new double[n];
  double *r = new double[n];
  double *s = new double[n];
  double *q = new double[n];

  // Stop condition limit.
  double tol_norm_b = tol * norm_b;

  // Setup the residual, r = b - A x
  A_residual(x, r);

  // Preconditioning, d = M^-1 r
  preconditioner_solve(r, d);

  // delta = r^T d
  double delta = ds_vector_dot(n, r, d);

  // residual norm
  double sum_r2 = 0.0;
  for (int i = 0; i < n; ++i) { sum_r2 += r[i] * r[i]; }
  double norm_r = sqrt(sum_r2);
  res.residual_norms.push_back(norm_r);

  if (verbose) {
    printf("[PCG] Solving %i x %i problem.\n", n, n);
    printf("    Goal |r| < (tol) |b| = %e,  delta = %e\n", tol_norm_b, delta);
  }

  int iter;
  for (iter = 1; iter < max_iter; ++iter) {
    bool full_residual = iter % step_r == 0;
    pcg_step(n, A_product, A_residual, preconditioner_solve,
        d, q, r, s, x, full_residual, &delta);

    // compute residual norm
    sum_r2 = 0.0;
    for (int i = 0; i < n; ++i) { sum_r2 += r[i] * r[i]; }
    norm_r = sqrt(sum_r2);
    res.residual_norms.push_back(norm_r);

    // Check stop condition.
    if (norm_r <= tol_norm_b) {
      if (verbose) {
          printf("    iter %i reached |r| = %e < tol |b| = %e\n", iter, norm_r, tol_norm_b);
      }
      break;
    }

    // Output progress.
    if (verbose) {
      printf("    iter %i, |r| %e, delta %e\n", iter, norm_r, delta);
    }
  }

  // don't forget to free the work arrays.
  delete [] d;
  delete [] r;
  delete [] q;
  delete [] s;

  res.success = true;

  // max iter warning
  if (iter == max_iter) {
    res.success = false;
    puts("[WARNING] PCG reached max iteration before stop condition.");
  }

  // Update res before returning
  res.num_iters = iter;
  return res;
}

#endif
