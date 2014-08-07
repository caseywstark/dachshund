/*
2014, the dachshund authors.
*/

#include <cmath>

#include "dachshund.h"
#include "test_utils.h"

const int num_pixels = 2000;
const std::string pixel_data_path = "tests/small_pixel_data.bin";

struct SmallWFX
{
    SignalCovarParams *s_params;

    Pixel *pixels;
    NPixel *npixels;
    WPixel *wpixels;

    WFX_SN *wfx_sn;
    WFX_W *wfx_w;

    SmallWFX() {
        s_params = new SignalCovarParams(0.05, 1.0, 1.0);

        // pixel setup
        pixels = new Pixel[num_pixels];
        read_pixel_data(pixel_data_path, num_pixels, pixels);

        npixels = new NPixel[num_pixels];
        pixel_to_npixel(num_pixels, pixels, npixels);

        wpixels = new WPixel[num_pixels];
        pixel_to_wpixel(num_pixels, pixels, wpixels);

        wfx_sn = new WFX_SN(num_pixels, npixels, s_params);
        wfx_w = new WFX_W(num_pixels, wpixels, s_params);
    }

    ~SmallWFX() {
        delete s_params;
        delete [] pixels;
        delete [] npixels;
        delete [] wpixels;
        delete wfx_sn;
        delete wfx_w;
    }
};

TEST_FIXTURE(SmallWFX, sn_vs_w)
{
    puts("[TEST] x solve: different x expressions.");

    double *x_sn = new double[num_pixels];
    double *x_w = new double[num_pixels];

    // false for verbose, true for exact
    wfx_sn->solve_cf(x_sn, false, true);
    wfx_w->solve_cf(x_w, false, true);

    // test for close
    for (int i = 0; i < num_pixels; ++i) {
        CHECK_DTOL(x_sn[i], x_w[i], 0.0, 1.0e-8);
    }

    delete [] x_sn;
    delete [] x_w;
}

TEST_FIXTURE(SmallWFX, Aij_lookup_s_vs_exact_s)
{
    puts("[TEST] x solve: A_ij with lookup S vs. exact S.");

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

    CHECK(xerr_rms < 1.0e-3 * x_rms);

}

TEST_FIXTURE(SmallWFX, lookup_s_vs_exact_s)
{
    puts("[TEST] x solve: exact solve with lookup S vs. exact S.");

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

    CHECK(xerr_rms < 1.0e-2 * x_rms);

    delete [] x_lu;
    delete [] x_exact;
}


TEST_FIXTURE(SmallWFX, solve_pcg_vs_cf)
{
    puts("[TEST] x solve: PCG with lookup S vs. CF with exact S.");
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

    CHECK(rms < 0.01);

    delete [] x_pcg;
    delete [] x_cf;
}

/*
TEST_FIXTURE(SmallWFX, pcg_sn_vs_w)
{
    double *x_chol = new double[num_pixels];
    double *x_sn = new double[num_pixels];
    for (int i = 0; i < num_pixels; ++i) { x_sn[i] = 0.0; }
    double *x_w = new double[num_pixels];
    for (int i = 0; i < num_pixels; ++i) { x_w[i] = 0.0; }

    // get exact solution
    wfx_sn->chol_solve(x_chol, false);

    // set pcg params
    pcg_params.max_iter = 1000;
    pcg_params.step_r = 50;
    pcg_params.tol = 1.0e-3;

    puts("\n[TEST] SN PCG.");
    wfx_sn->pcg_solve(x_sn, &pcg_params, true);

    puts("\n[TEST] W PCG.");
    wfx_w->pcg_solve(x_w, &pcg_params, true);

    puts("\n[TEST] SN vs. W PCG solutions");
    for (int i = 0; i < 50; ++i) {
        printf("i %i \t PCG %e \t CHOL %e \t error %f\n", i, x_pcg[i], x_chol[i], x_pcg[i] / x_chol[i] - 1.0);
    }
}
*/
