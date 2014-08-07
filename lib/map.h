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

//
// signal covar
//

const double default_gaussian_table_x_max = 4.0;
const int default_gaussian_table_n = 200;

class SignalCovarParams {
  public:
    double var_s;
    double gaussian_const_perp;
    double gaussian_const_para;
    int func_table_n;
    double func_table_dx_inv;
    double *func_table;

    SignalCovarParams(const double var_s, const double l_perp,
        const double l_para)
        : var_s(var_s)
    {
        gaussian_const_perp = 1.0 / (2.0 * l_perp * l_perp);
        gaussian_const_para = 1.0 / (2.0 * l_para * l_para);

        // enforced defaults
        const double fx_max = default_gaussian_table_x_max
            * default_gaussian_table_x_max / 2.0;
        func_table_n = default_gaussian_table_n;
        const double ft_dx = fx_max / (func_table_n - 1);

        func_table_dx_inv = 1.0 / ft_dx;
        func_table = new double[func_table_n];

        // Iterate over points x_i = dx * i, filling gaussian function table.
        // Note that we use exp(-x) here, since we want to pass x^2
        // when we use the table.
        for (int i = 0; i < func_table_n; ++i) {
            double xi = ft_dx * i;
            func_table[i] = exp(-xi);
        }
    }

    ~SignalCovarParams()
    {
        delete [] func_table;
    }
};

/*
Linear interpolation, given x and table data.
We assume x >= 0 and the table data is f_i = f(x_i) = f(dx * i).
If x is larger than the last tabulated point, we return 0.
*/
static inline double
linterp(const double x, const int n, const double dx_inv,
    const double * const f)
{
    // in index space.
    double xx = x * dx_inv;
    // decimate to get low edge.
    int ix = xx;
    // Check if we are outside table range.
    if (ix > n - 2) {
        return 0.0;
    }
    // the linterp expression.
    return f[ix] * (ix + 1 - xx) + f[ix+1] * (xx - ix);
}

static inline double
signal_covar(const double x_perp_2, const double x_para_2,
    const SignalCovarParams * const p)
{
    const double ga_perp = p->gaussian_const_perp * x_perp_2;
    const double ga_para = p->gaussian_const_para * x_para_2;

    const double g_perp = linterp(ga_perp, p->func_table_n,
        p->func_table_dx_inv, p->func_table);
    const double g_para = linterp(ga_para, p->func_table_n,
        p->func_table_dx_inv, p->func_table);

    const double s = p->var_s * g_perp * g_para;
    return s;
}

static inline double
signal_covar_exact(const double x_perp_2, const double x_para_2,
    const SignalCovarParams * const p)
{
    const double ga_perp = p->gaussian_const_perp * x_perp_2;
    const double ga_para = p->gaussian_const_para * x_para_2;

    double s = p->var_s * exp(-ga_perp) * exp(-ga_para);
    return s;
}

// PCG data structs

struct pcg_params_s {
    int max_iter;
    int step_r;
    double tol;
};

typedef struct pcg_params_s PCGParams;

struct pcg_result_s {
    int num_iters;
    std::vector<double> residual_norms;
};

typedef struct pcg_result_s PCGResult;

//
// Wiener Filtering vector versions.
//

// The default vector form is x = A^{-1} b = (S + N)^{-1} d
class WFX_SN {
  public:
    const int num_pixels;
    const NPixel * const pixels;
    const SignalCovarParams * const s_params;

    WFX_SN(const int num_pixels, const NPixel * const pixels,
        const SignalCovarParams * const s_params)
        : num_pixels(num_pixels), pixels(pixels), s_params(s_params)
    { }

    double A(const int i, const int j, const bool exact=false);

    // x = A~^{-1} b where A~ is the preconditioner matrix
    void preconditioner_solve(const double * const b, double * const x);
    // b = A a
    void A_product(const double * const a, double * const b);
    // r = b - A x
    void A_residual(const double * const x, double * const r);
    // Just |b|
    double norm_b();

    // One iteration of the PCG.
    // d, q, and s are misc. work arrays.
    // r is the residual. x is the current solution.
    // full_residual controls updating the residual with the approximate or
    // full expression.
    // delta is the step size which needs to be used for the next iteration.
    void pcg_step(double * const d, double * const q, double * const r,
        double * const s, double * const x, const bool full_residual,
        double * delta);

    // Runs PCG.
    PCGResult solve_pcg(double * const x, const PCGParams * const pcg_params,
        const bool verbose=true);

    // Runs a direct cholesky factorization.
    void solve_cf(double * const x, const bool verbose=true, const bool exact=false);
};

class WFX_W {
  public:
    const int num_pixels;
    const WPixel * const pixels;
    const SignalCovarParams * const s_params;

    WFX_W(const int num_pixels, const WPixel * const pixels,
        const SignalCovarParams * const s_params)
        : num_pixels(num_pixels), pixels(pixels), s_params(s_params)
    { }

    double A(const int i, const int j, bool exact=false);

    // x = A~^{-1} b where A~ is the preconditioner matrix
    void preconditioner_solve(const double * const b, double * const x);
    // b = A a
    void A_product(const double * const a, double * const b);
    // r = b - A x
    void A_residual(const double * const x, double * const r);
    // Just |b|
    double norm_b();

    // One iteration of the PCG.
    // d, q, and s are misc. work arrays.
    // r is the residual. x is the current solution.
    // full_residual controls updating the residual with the approximate or
    // full expression.
    // delta is the step size which needs to be used for the next iteration.
    void pcg_step(double * const d, double * const q, double * const r,
        double * const s, double * const x, const bool full_residual,
        double * delta);

    // Runs PCG.
    PCGResult solve_pcg(double * const x, const PCGParams * const pcg_params,
        const bool verbose=true);

    // Runs a direct cholesky factorization.
    void solve_cf(double * const x, const bool verbose=true, const bool exact=false);
};

void
smp_product(const int num_pixel_points, const Point * const pixel_coords,
        const int num_map_points, const Point * const map_coords,
        const SignalCovarParams * const s_params,
        const double * const x, double * const m);

/*
class Reconstruction {
  public:
};
*/

#endif
