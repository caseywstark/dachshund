/*
2015, the dachshund authors.
*/

#include <cmath>

#include "gtest/gtest.h"

#include "dachshund.h"
#include "test_utils.h"

class MapCovarTest : public testing::Test {
  protected:
    const static int num_pixels = 2000;

    const double lx = 10.0;
    const double ly = 10.0;
    const double lz = 10.0;
    const int map_nx = 1;
    const int map_ny = 1;
    const int map_nz = 5;
    const int num_map_points = 5;

    virtual void SetUp() {
        const std::string pixel_data_path = "tests/small_pixel_data.bin";

        s_params = new SignalCovarParams(0.05, 2.0, 2.0);
        pcg_params = {num_pixels, 1000, 50, 1e-2, false};

        pixels = new Pixel[num_pixels];
        read_pixel_data(pixel_data_path, num_pixels, pixels);
        npixels = new NPixel[num_pixels];
        pixel_to_npixel(num_pixels, pixels, npixels);

        map_coords = new Point[num_map_points];
        init_uniform_map_coords(lx, ly, lz, map_nx, map_ny, map_nz, map_coords);

        map_covar_diag = new MapCovarDiag(num_pixels, npixels, num_map_points,
            map_coords, s_params);

        sigma_m = new double[num_pixels];
    }

    virtual void TearDown() {
        delete s_params;
        delete [] pixels;
        delete [] npixels;
        delete [] map_coords;
        delete [] sigma_m;
        delete map_covar_diag;
    }

    SignalCovarParams* s_params;
    Pixel* pixels;
    NPixel* npixels;
    Point* map_coords;
    double* sigma_m;
    PCGParams pcg_params;
    MapCovarDiag* map_covar_diag;

};

TEST_F(MapCovarTest, map_covar_solve) {
    map_covar_diag->solve_pcg(pcg_params, sigma_m);

    FILE *f = fopen("map_covar.bin", "w");
    fwrite(sigma_m, sizeof(double), num_map_points, f);
    fclose(f);
}

TEST_F(MapCovarTest, noise_dom_solve) {
    sigma_m_noise_dom(num_pixels, npixels, num_map_points, map_coords,
        s_params, sigma_m);
}
