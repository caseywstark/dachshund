/*
2014, the dachshund authors.
*/

#include <cmath>

#include "dachshund.h"
#include "test_utils.h"

const int mock_n = 1718;
const std::string mock_data_path = "tests/pilot_mock_data.bin";

struct Mock
{
    SignalCovarParams *s_params;

    Pixel *mock_pixels;
    WPixel *mock_wpixels;
    WFX_W *wfx;

    Mock() {
        // signal covar
        s_params = new SignalCovarParams(0.10, 3.0, 2.0);

        // pixel setup
        mock_pixels = new Pixel[mock_n];
        read_pixel_data(mock_data_path, mock_n, mock_pixels);

        // wfx setup
        mock_wpixels = (WPixel *) mock_pixels;
        pixel_to_wpixel(mock_n, mock_pixels, mock_wpixels);
        wfx = new WFX_W(mock_n, mock_wpixels, s_params);
    }

    ~Mock() {
        delete s_params;
        delete [] mock_pixels;
        delete wfx;
    }
};


TEST_FIXTURE(Mock, mock_pcg_vs_cf)
{
    puts("[TEST] Pilot mock map solution PCG vs. CF");

    double *x_cf = new double[mock_n];
    wfx->solve_cf(x_cf, false, true);

    double *x_pcg = new double[mock_n];
    for (int i = 0; i < mock_n; ++i) { x_pcg[i] = 0.0; }
    wfx->solve_cf(x_pcg, false, true);

    // setup map coordinates
    const double lx = 6.0;
    const double ly = 14.0;
    const double lz = 80.0;
    const int map_nx = lx / 0.5;
    const int map_ny = ly / 0.5;
    const int map_nz = lz / 0.5;
    const int num_map_points = map_nx * map_ny * map_nz;
    Point * const map_coords = new Point[num_map_points];
    init_uniform_map_coords(lx, ly, lz, map_nx, map_ny, map_nz, map_coords);

    // pixel points
    const Point * const pixel_coords =
        points_from_pixels_alloc(mock_n, mock_pixels);

    // map multiplies
    double *m_cf = new double[num_map_points];
    double *m_pcg = new double[num_map_points];

    smp_product(mock_n, pixel_coords, num_map_points, map_coords,
        s_params, x_cf, m_cf);
    smp_product(mock_n, pixel_coords, num_map_points, map_coords,
        s_params, x_pcg, m_pcg);

    // compute map difference rms
    double md2_sum = 0.0;
    for (int i = 0; i < num_map_points; ++i) {
        double m_diff = m_pcg[i] / m_cf[i] - 1.0;
        md2_sum += m_diff * m_diff;
    }

    double md2_mean = md2_sum / num_map_points;
    double md_rms = sqrt(md2_mean);

    CHECK(md_rms < 0.01);
}


/*
TEST_FIXTURE(Mock, zero_weights)
{
    // Exact cholesky solve with W version and zeros.
    // then compare to solutions with N version and different values for
    // the masked large noise.
    WPixel *wpixels = new WPixel[num_pixels];
    pixel_to_wpixel(num_pixels, pixels, wpixels);

    // exact solution with masking
    WFX_W *wfx_w = new WFX_W(num_pixels, wpixels, s_params);
    double *x_mask = new double[num_pixels];
    wfx_w->chol_exact_solve(x_mask);

    // Create noise version
    NPixel *npixels = new NPixel[num_pixels];

    double n_max = 100.0;
    for (int i = 0; i < num_pixels; ++i) {
        if (pixels[i].n > n_max) {
            pixels[i].n = n_max;
        }
    }

    pixel_to_npixel(num_pixels, pixels, npixels);

    WFX_SN *wfx_sn = new WFX_SN(num_pixels, npixels, s_params);
    double *x_n1 = new double[num_pixels];
    wfx_sn->chol_exact_solve(x_n1);

    for (int i = 0; i < num_pixels; ++i) {
        printf("i %i \t MASK %e \t N1 %e \t error %f\n", i, x_mask[i], x_n1[i], x_n1[i] / x_mask[i] - 1.0);
    }
}
*/
