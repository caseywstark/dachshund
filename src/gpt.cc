/*
2014, the dachshund authors.
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "gpt.h"

/*
Bilinear interpolation, given (x, y) and 2d table data.

We assume x, y >= 0. If x or y are larger than the last tabulated point, we
return 0.

Parameters
x, y : the eval points
nx, ny : the number of samples
dx, dy : the sample point separation
f_table : array containing the function value at sample points.
*/
double
bilinterp(const double x, const double y, const int nx, const int ny,
    const double dx, const double dy, const double * const f_table)
{
    // x and y in index space.
    double xx = x / dx;
    double yy = y / dy;

    // the lower point index.
    int ix = xx;
    int iy = yy;

    // Check if we are outside table range.
    if (ix > nx - 1 || iy > ny - 1) {
        return 0;
    }

    // The weights.
    double wx0 = (ix+1) - xx;
    double wx1 = xx - ix;
    double wy0 = (iy+1) - yy;
    double wy1 = yy - iy;

    // Grab corner values.
    double f00 = f_table[ix * ny + iy];
    double f01 = f_table[ix * ny + iy+1];
    double f10 = f_table[(ix+1) * ny + iy];
    double f11 = f_table[(ix+1) * ny + iy+1];

    // bilinear expression.
    double f = f00 * wx0 * wy0 + f01 * wx0 * wy1 + f10 * wx1 * wy0
        + f11 * wx1 * wy1;

    // DEBUG
    /*
    printf("x y %g %g   xx yy %g %g   ix iy %i %i   w %g %g %g %g\n",
        x, y, xx, yy, ix, iy, wx0, wx1, wy0, wy1);
    */

    return f;
}

GPT::GPT(const int n, const double dx)
{
    this->n = n;
    this->dx = dx;

    // Allocate table
    int nn = n * n;
    table = new double[nn];

    // Set table values.
    int ij = 0;
    for (int i = 0; i < n; ++i) {
        double xi = i * dx;
        double g1 = exp(-xi * xi);
        for (int j = 0; j < n; ++j) {
            double xj = j * dx;
            double g2 = exp(-xj * xj);
            // set t[ij] and *then* increment ij.
            table[ij++] = g1 * g2;
        }
    }
}

GPT::~GPT()
{
    delete table;
}

double
GPT::operator()(const double x, const double y)
{
    double f = bilinterp(x, y, n, n, dx, dx, table);
    return f;
}
