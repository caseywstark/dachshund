/*
2014, the dachshund authors.
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "util.h"

void
read_pixel_data(const std::string path, const int num_pixels,
    Pixel * const pixels)
{
    FILE *infile = NULL;
    infile = fopen(path.c_str(), "r");
    if (infile == NULL) {
        printf("Could not open file %s\n", path.c_str());
        exit(1);
    }
    fread(pixels, sizeof(Pixel), num_pixels, infile);
    fclose(infile);
}

void
init_uniform_map_coords(const double lx, const double ly, const double lz,
    const int nx, const int ny, const int nz, Point * const coords)
{
    const double dx = lx / nx;
    const double dy = ly / ny;
    const double dz = lz / nz;

    for (int ix = 0; ix < nx; ++ix) {
        const double x = dx * (ix + 0.5);
        for (int iy = 0; iy < ny; ++iy) {
            const double y = dy * (iy + 0.5);
            for (int iz = 0; iz < nz; ++iz) {
                const double z = dz * (iz + 0.5);
                int i = (ix * ny + iy) * nz + iz;
                coords[i].x = x;
                coords[i].y = y;
                coords[i].z = z;
            }
        }
    }
}

void
check_finite(const int n, const double * const v)
{
    for (int i = 0; i < n; ++i) {
        bool f = std::isfinite(v[i]);
        if (!f) {
            printf("non-finite value v[%i] = %e\n", i, v[i]);
        }
    }
}
