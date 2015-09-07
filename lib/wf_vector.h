//
// Wiener Filter vector versions (the x in x = A^{-1} b).
//

#ifndef __DS_WF_VECTOR_H__
#define __DS_WF_VECTOR_H__

#include <functional>

#include "pcg.h"
#include "pixel.h"
#include "signal_covar.h"

class WFVector {
  public:
    WFVector(const int num_pixels, const SignalCovarParams* const s_params)
        : n(num_pixels), s_params(s_params) { }

    const int n;
    const SignalCovarParams* const s_params;

    double A(const int i, const int j, const bool exact);
    void preconditioner_solve(const double* const b, double * const x);
    void A_product(const double* const a, double* const b);
    void A_residual(const double* const x, double* const r);
    void norm_b();

    PCGResult solve_pcg(const PCGParams& pcg_params, double* const x);
    void solve_cf(double* const x, const bool verbose, const bool exact);
};

// The default vector form, x = (S + N)^{-1} d
class WFVectorSN : public WFVector {
  public:
    const NPixel* const pixels;

    WFVectorSN(const int num_pixels, const NPixel* const pixels,
        const SignalCovarParams* const s_params)
        : WFVector(num_pixels, s_params), pixels(pixels) {
      this_pc_solve = std::bind(&WFVectorSN::preconditioner_solve, this,
          std::placeholders::_1, std::placeholders::_2);
      this_A_product = std::bind(&WFVectorSN::A_product, this,
          std::placeholders::_1, std::placeholders::_2);
      this_A_residual = std::bind(&WFVectorSN::A_residual, this,
          std::placeholders::_1, std::placeholders::_2);
    }

    double A(const int i, const int j, const bool exact);
    double norm_b();

    void preconditioner_solve(const double* const b, double* const x);
    void A_product(const double* const a, double* const b);
    void A_residual(const double* const x, double* const r);

    std::function<void(const double* const, double* const)> this_pc_solve;
    std::function<void(const double* const, double* const)> this_A_product;
    std::function<void(const double* const, double* const)> this_A_residual;

    PCGResult solve_pcg(const PCGParams& pcg_params, double* const x);
    void solve_cf(double* const x, const bool verbose, const bool exact);
};

//
// WF version x = w (w S w + I)^{-1} w d
//   A_ij = w_i w_j S_ij + delta_ij
//   b_i = w_i d_i
//
class WFVectorW : public WFVector {
  public:
    const WPixel* const pixels;

    WFVectorW(const int num_pixels, const WPixel* const pixels,
        const SignalCovarParams* const s_params)
        : WFVector(num_pixels, s_params), pixels(pixels) {
      this_pc_solve = std::bind(&WFVectorW::preconditioner_solve, this,
          std::placeholders::_1, std::placeholders::_2);
      this_A_product = std::bind(&WFVectorW::A_product, this,
          std::placeholders::_1, std::placeholders::_2);
      this_A_residual = std::bind(&WFVectorW::A_residual, this,
          std::placeholders::_1, std::placeholders::_2);
    }

    double A(const int i, const int j, const bool exact);
    double norm_b();

    void preconditioner_solve(const double* const b, double* const x);
    void A_product(const double* const a, double* const b);
    void A_residual(const double* const x, double* const r);

    std::function<void(const double* const, double* const)> this_pc_solve;
    std::function<void(const double* const, double* const)> this_A_product;
    std::function<void(const double* const, double* const)> this_A_residual;

    PCGResult solve_pcg(const PCGParams& pcg_params, double* const x);
    void solve_cf(double* const x, const bool verbose, const bool exact);
};

#endif
