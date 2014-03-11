/*
2014, the dachshund authors.
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "gt.h"

/*
Linear interpolation, given x and table data.

We assume x >= 0. If x is larger than the last tabulated point, we return 0.

Parameters
x : the eval point
nx : the number of samples
dx: the sample point separation
f_table : array containing the function value at sample points.
*/
double
linterp(const double x, const int nx, const double dx,
    const double * const f_table)
{
    // index space.
    double xx = x / dx;
    // the lower point index.
    int ix = xx;
    // Check if we are outside table range.
    if (ix > nx - 2) {
        return 0.0;
    }

    // The weights.
    double w0 = (ix+1) - xx;
    double w1 = xx - ix;

    // The bordering function values.
    double f0 = f_table[ix];
    double f1 = f_table[ix+1];

    double f = f0 * w0 + f1 * w1;

    return f;
}

GT::GT(const int n, const double dx)
{
    this->n = n;
    this->dx = dx;

    // Allocate table
    table = new double[n];

    // Set table values.
    for (int i = 0; i < n; ++i) {
        double xi = i * dx;
        table[i] = exp(-xi * xi);
    }
}

GT::~GT()
{
    delete [] table;
}

double
GT::f(const double x)
{
    double f = linterp(x, n, dx, table);
    return f;
}
