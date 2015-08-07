/*
2014, the dachshund authors.
*/

#include <cmath>

#include "gtest/gtest.h"

#include "dachshund.h"
#include "test_utils.h"

class SolveRegressionTest : public testing::Test {
  protected:
    const static int num_pixels = 2000;

    virtual void SetUp() {
        const std::string pixel_data_path = "tests/small_pixel_data.bin";

        s_params = new SignalCovarParams(0.05, 1.0, 1.0);

        pixels = new Pixel[num_pixels];
        read_pixel_data(pixel_data_path, num_pixels, pixels);
        npixels = new NPixel[num_pixels];
        pixel_to_npixel(num_pixels, pixels, npixels);
        wpixels = new WPixel[num_pixels];
        pixel_to_wpixel(num_pixels, pixels, wpixels);

        wfx_sn = new WFX_SN(num_pixels, npixels, s_params);
        wfx_w = new WFX_W(num_pixels, wpixels, s_params);
    }

    virtual void TearDown() {
        delete s_params;
        delete [] pixels;
        delete [] npixels;
        delete [] wpixels;
        delete wfx_sn;
        delete wfx_w;
    }

    SignalCovarParams *s_params;
    Pixel *pixels;
    NPixel *npixels;
    WPixel *wpixels;
    WFX_SN *wfx_sn;
    WFX_W *wfx_w;
};

TEST_F(SolveRegressionTest, different_x_expressions) {
    double *x_sn = new double[num_pixels];
    double *x_w = new double[num_pixels];

    // false for verbose, true for exact
    wfx_sn->solve_cf(x_sn, false, true);
    wfx_w->solve_cf(x_w, false, true);

    // test for close
    for (int i = 0; i < num_pixels; ++i) {
        EXPECT_EQTOL(x_sn[i], x_w[i], 0.0, 1.0e-8);
    }

    delete [] x_sn;
    delete [] x_w;
}

// This isn't really testing the solve, but it's simpler to use the problem
// setup to check A_ij than to separate it.
TEST_F(SolveRegressionTest, A_using_lookup_S_error_rms) {
    // compute A_ij error rms
    double x2_sum = 0.0;
    double xerr2_sum = 0.0;
    for (int i = 0; i < num_pixels; ++i) {
        for (int j = 0; j < num_pixels; ++j) {
            double aij_exact = wfx_sn->A(i, j, true);
            double aij_lu = wfx_sn->A(i, j, false);

            double x = aij_exact;
            double xerr = aij_lu - aij_exact;
            x2_sum += x * x;
            xerr2_sum += xerr * xerr;
        }
    }

    const int nn = num_pixels * num_pixels;
    double x2_mean = x2_sum / nn;
    double xerr2_mean = xerr2_sum / nn;
    double x_rms = sqrt(x2_mean);
    double xerr_rms = sqrt(xerr2_mean);

    EXPECT_LT(xerr_rms, 1.0e-3 * x_rms);
}

TEST_F(SolveRegressionTest, solve_using_lookup_S_vs_exact_S) {
    double *x_lu = new double[num_pixels];
    double *x_exact = new double[num_pixels];

    wfx_sn->solve_cf(x_lu, false, false);
    wfx_sn->solve_cf(x_exact, false, true);

    // compute x error rms
    double x2_sum = 0.0;
    double xerr2_sum = 0.0;
    for (int i = 0; i < num_pixels; ++i) {
        double x = x_exact[i];
        double xerr = x_lu[i] - x_exact[i];
        x2_sum += x * x;
        xerr2_sum += xerr * xerr;
    }

    const int nn = num_pixels;
    double x2_mean = x2_sum / nn;
    double xerr2_mean = xerr2_sum / nn;
    double x_rms = sqrt(x2_mean);
    double xerr_rms = sqrt(xerr2_mean);

    EXPECT_LT(xerr_rms, 1.0e-2 * x_rms);

    delete [] x_lu;
    delete [] x_exact;
}

// This is like the full error case.
TEST_F(SolveRegressionTest, PCG_lookup_S_vs_CF_exact_S) {
    PCGParams pcg_params;
    pcg_params.max_iter = 1000;
    pcg_params.step_r = 10;
    pcg_params.tol = 1.0e-4;

    double *x_pcg = new double[num_pixels];
    for (int i = 0; i < num_pixels; ++i) { x_pcg[i] = 0.0; }
    double *x_cf = new double[num_pixels];

    wfx_w->solve_cf(x_cf, false, false);
    wfx_w->solve_pcg(x_pcg, &pcg_params, false);

    // compute x error rms
    double x2_sum = 0.0;
    for (int i = 0; i < num_pixels; ++i) {
        double x_err = x_pcg[i] / x_cf[i] - 1.0;
        x2_sum += x_err * x_err;
    }
    double x2_mean = x2_sum / num_pixels;
    double rms = sqrt(x2_mean);

    EXPECT_LT(rms, 0.01);

    delete [] x_pcg;
    delete [] x_cf;
}
