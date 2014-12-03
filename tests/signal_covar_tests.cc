/*
2014, the dachshund authors.
*/

#include <cmath>

#include "catch.hpp"

#include "dachshund.h"
#include "test_utils.h"

TEST_CASE("S lookup vs. exact 1D", "[signal-covar]")
{
    const double s_tol = 0.001;
    const double l = 0.5;
    const double x_max = default_gaussian_table_x_max * l;
    const double dx_samp = 0.03 * l;

    SignalCovarParams *p = new SignalCovarParams(1.0, l, l);

    double x = 0.0;
    while (x < x_max) {
        // eval fast and exact S
        double sf = signal_covar(x*x, 0.0, p);
        double se = exp(- (x * x) / 2.0 / (l * l) );
        //printf("signal covar: x %e   Sf %e   Se %e\n", x, sf, se);
        REQUIRE( dtol(sf, se, 0.0, s_tol) );
        x += dx_samp;
    }
}

TEST_CASE("S lookup vs. exact", "[signal-covar]")
{
    const double s_tol = 0.001;

    const double l0 = 0.8;
    const double x0_max = default_gaussian_table_x_max * l0;
    const double l1 = 1.2;
    const double x1_max = default_gaussian_table_x_max * l1;

    const double dx0 = 0.04 * l0;
    const double dx1 = 0.04 * l1;

    SignalCovarParams *p = new SignalCovarParams(1.0, l0, l1);

    const int n0 = x0_max / dx0;
    const int n1 = x1_max / dx1;

    for (int i0 = 0; i0 < n0; ++i0) {
        double x0 = i0 * dx0;
        for (int i1 = 0; i1 < n1; ++i1) {
            double x1 = i1 * dx1;

            // eval fast and exact S
            double sf = signal_covar(x0*x0, x1*x1, p);
            if (sf == 0.0) { continue; } // might hit some of these points.

            double se = exp(-(x0 * x0) / 2.0 / (l0 * l0) -(x1 * x1) / 2.0 / (l1 * l1));

            //printf("signal covar: x0 %e x1 %e   Sf %e   Se %e\n", x0, x1, sf, se);
            REQUIRE( dtol(sf, se, 0.0, s_tol) );
        }
    }
}


TEST_CASE("S lookup table cutoff", "[signal-covar]")
{
    const double l = rng();
    SignalCovarParams *p = new SignalCovarParams(1.0, l, l);

    // we usually want s to cut out at x >= 3 sigma.
    double x, sf;

    x = (default_gaussian_table_x_max - 0.01) * l;
    sf = signal_covar(x*x, 0.0, p);
    REQUIRE(sf > 0.0);

    x = (default_gaussian_table_x_max + 0.01) * l;
    sf = signal_covar(x*x, 0.0, p);
    REQUIRE(sf == 0.0);
}
