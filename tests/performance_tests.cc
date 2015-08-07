/*
2015, the dachshund authors.
*/

#include <cmath>

#include "gtest/gtest.h"

#include "dachshund.h"
#include "test_utils.h"

class PerformanceRegressionTest : public testing::Test {
  protected:
    const static int num_pixels = 2000;

    virtual void SetUp() {
        const std::string pixel_data_path = "tests/small_pixel_data.bin";

        s_params = new SignalCovarParams(0.05, 2.0, 2.0);

        pixels = new Pixel[num_pixels];
        read_pixel_data(pixel_data_path, num_pixels, pixels);
        npixels = new NPixel[num_pixels];
        pixel_to_npixel(num_pixels, pixels, npixels);
        wpixels = new WPixel[num_pixels];
        pixel_to_wpixel(num_pixels, pixels, wpixels);

        wfx_sn = new WFX_SN(num_pixels, npixels, s_params);
        wfx_w = new WFX_W(num_pixels, wpixels, s_params);

        x = new double[num_pixels];
    }

    virtual void TearDown() {
        delete s_params;
        delete [] pixels;
        delete [] npixels;
        delete [] wpixels;
        delete wfx_sn;
        delete wfx_w;
    }

    Timer timer;
    SignalCovarParams *s_params;
    Pixel *pixels;
    NPixel *npixels;
    WPixel *wpixels;
    WFX_SN *wfx_sn;
    WFX_W *wfx_w;
    double *x;
};

TEST_F(PerformanceRegressionTest, wfx_sn_solve_time) {
    timer.reset();
    wfx_sn->solve_cf(x, false, false);
    double elapsed = timer.elapsed();
    printf("    Took %.2f s (typically 0.7 s on 2.2 GHz)\n", elapsed);
    EXPECT_LT(elapsed, 1.0);
}

TEST_F(PerformanceRegressionTest, wfx_w_solve_time) {
    timer.reset();
    wfx_w->solve_cf(x, false, false);
    double elapsed = timer.elapsed();
    printf("    Took %.2f s (typically 0.7 s on 2.2 GHz)\n", elapsed);
    EXPECT_LT(elapsed, 1.0);
}

TEST_F(PerformanceRegressionTest, map_multiply_time) {
    // This is a waste, but small enough.
    wfx_sn->solve_cf(x, false, false);

    const Point * const pixel_coords =
        points_from_pixels_alloc(num_pixels, pixels);

    const int n_map = 50;
    const int num_map_points = n_map * n_map * n_map;
    Point * const map_coords = new Point[num_map_points];
    double *m = new double[num_map_points];

    init_uniform_map_coords(1.0, 1.0, 1.0, n_map, n_map, n_map, map_coords);

    timer.reset();
    smp_product(num_pixels, pixel_coords, num_map_points, map_coords,
        s_params, x, m);
    double elapsed = timer.elapsed();
    printf("    Took %.2f s (typically 0.9 s on 2.2 GHz)\n", elapsed);
    EXPECT_LT(elapsed, 2.0);
}

// TODO: Add one for covariance calculation.
