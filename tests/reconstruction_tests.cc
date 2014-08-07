/*
2014, the dachshund authors.
*/

#include "dachshund.h"
#include "test_utils.h"

/*
struct RandomPixels
{
    RandomPixels() {
        const double l = 1.0;
        const int n = 8;

        num_pixels = n*n*n;
        num_map_points = num_pixels;

        pixels = new Pixel[num_pixels];
        map_coords = new Point[num_map_points];

        for (int ix = 0; ix < n; ++ix) {
            for (int iy = 0; iy < n; ++iy) {
                for (int iz = 0; iz < n; ++iz) {
                    int i = (ix * n + iy) * n + iz;
                    pixels[i].x = l / n * (ix + 0.5);
                    pixels[i].y = l / n * (iy + 0.5);
                    pixels[i].z = l / n * (iz + 0.5);
                    double n = 1.0e-3;
                    pixels[i].N = n*n;
                    pixels[i].d = rng();
                }
            }
        }

        init_uniform_map_coords(l, l, l, n, n, n, map_coords);

        signal_params = new FastSignalParams(1.0, 1.0e-3, 1.0e-3);

        pcg_params.max_iter = 100;
        pcg_params.step_r = 50;
        pcg_params.tol = 1.0e-6;

        map = new Map(num_pixels, pixels, num_map_points, map_coords,
            signal_params, &pcg_params);

        m = new double[num_map_points];
    }

    ~RandomPixels() {
        delete [] map_coords;
        delete [] pixels;
        delete map;
        delete [] m;
    }

    int num_pixels;
    Pixel *pixels;
    int num_map_points;
    Point *map_coords;
    FastSignalParams *signal_params;
    PCGParams pcg_params;
    Map *map;
    double *m;
};

TEST_FIXTURE(RandomPixels, NoExtent)
{
    map->solve(m, false);

    for (int i = 0; i < num_pixels; ++i) {
        double d = pixels[i].d;
        //printf("i %i   d %e   m %e\n", i, d, m[i]);
        CHECK_DTOL(m[i], d, 0.0, 1.0e-6);
    }
}

*/
