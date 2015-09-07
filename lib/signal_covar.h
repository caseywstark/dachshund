
#ifndef __DS_SIGNAL_COVAR_H__
#define __DS_SIGNAL_COVAR_H__

#include <cmath>

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
      : var_s(var_s) {
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

  ~SignalCovarParams() {
    delete [] func_table;
  }
};

/*
Linear interpolation, given x and table data.
We assume x >= 0 and the table data is f_i = f(x_i) = f(dx * i).
If x is larger than the last tabulated point, we return 0.
*/
static inline double linterp(const double x, const int n, const double dx_inv,
    const double* const f) {
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

static inline double signal_covar(const double x_perp_2, const double x_para_2,
    const SignalCovarParams* const p) {
  const double ga_perp = p->gaussian_const_perp * x_perp_2;
  const double ga_para = p->gaussian_const_para * x_para_2;

  const double gp = linterp(ga_perp + ga_para, p->func_table_n,
      p->func_table_dx_inv, p->func_table);

  const double s = p->var_s * gp;
  return s;
}

static inline double signal_covar_exact(const double x_perp_2,
    const double x_para_2, const SignalCovarParams* const p) {
  const double ga_perp = p->gaussian_const_perp * x_perp_2;
  const double ga_para = p->gaussian_const_para * x_para_2;

  double s = p->var_s * exp( -(ga_perp + ga_para) );
  return s;
}

#endif
