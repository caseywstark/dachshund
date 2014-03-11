/*
2014, the dachshund authors.
*/

#include <math.h>
#include <cstdio>

int
main()
{
    int n = 4;
    int num_skewers = n*n;
    int num_pixels = n;
    int nn = num_skewers * num_pixels;

    double l = 1.0;
    double dx = l / n;

    // set coords
    double sx[num_skewers], sy[num_skewers];
    for (int ix = 0; ix < n; ++ix) {
        for (int iy = 0; iy < n; ++iy) {
            int i = ix * n + iy;
            sx[i] = dx * (ix + 0.5);
            sy[i] = dx * (iy + 0.5);
        }
    }

    double xc = 0.5;

    // set data and weights.
    double data[nn], w[nn];
    for (int ix = 0; ix < n; ++ix) {
        for (int iy = 0; iy < n; ++iy) {
            for (int iz = 0; iz < n; ++iz) {
                int i = (ix * n + iy) * n + iz;
                double x = dx * (ix + 0.5);
                double y = dx * (iy + 0.5);
                double z = dx * (iz + 0.5);

                double dx = x - xc;
                double dy = y - xc;
                double dz = z - xc;
                double xx = sqrt(dx*dx + dy*dy + dz*dz) / 0.3;

                data[i] = 1.0 * exp(-xx*xx);
                w[i] = 1.0;
            }
        }
    }

    FILE *skewer_file = fopen("skewer_x.bin", "w");
    fwrite(sx, sizeof(double), num_skewers, skewer_file);
    fclose(skewer_file);

    skewer_file = fopen("skewer_y.bin", "w");
    fwrite(sy, sizeof(double), num_skewers, skewer_file);
    fclose(skewer_file);

    FILE *data_file = fopen("pixel_data.bin", "w");
    fwrite(data, sizeof(double), nn, data_file);
    fclose(data_file);

    FILE *w_file = fopen("pixel_weights.bin", "w");
    fwrite(w, sizeof(double), nn, w_file);
    fclose(w_file);

    return 0;
}
