/*
2014, the dachshund authors.
*/

#include <cstdio>
#include <math.h>

#include "params.h"

#include "lum.h"

double
smd_element_func(const int i, const int j)
{
    // grid flat index i -> grid ix, iy, iz
    // probably a faster way to do this...
    int iz = i % p.nz;
    int iy = (i - iz) / p.nz % p.ny;
    int ix = (i - iz - iy * p.nz) / (p.ny * p.nz) % p.nx;

    double x1 = p.dx * (ix + 0.5);
    double y1 = p.dy * (iy + 0.5);
    double z1 = p.dz * (iz + 0.5);

    int isk = j / p.num_pixels;
    int iz2 = j % p.num_pixels;
    double x2 = p.skewer_x[isk];
    double y2 = p.skewer_y[isk];
    double z2 = p.dz_pix * (iz2 + 0.5);

    // Separation vector components.
    double dx12 = x2 - x1;
    double dy12 = y2 - y1;
    double dz12 = z2 - z1;

    // The gaussian terms.
    double x_perp = sqrt(dx12*dx12 + dy12*dy12) / p.l_perp;
    double x_para = fabs(dz12) / p.l_para;

    // Look up gaussian product.
    double gp = p.gt->f(x_perp) * p.gt->f(x_para);
    // Don't forget the variance.
    double sn = p.sigma * gp;

    return sn;
}


double
sddn_element_func(const int i, const int j)
{
    int isk, iz;
    isk = i / p.num_pixels;
    iz = i % p.num_pixels;
    double x1 = p.skewer_x[isk];
    double y1 = p.skewer_y[isk];
    double z1 = p.dz_pix * (iz + 0.5);

    isk = j / p.num_pixels;
    iz = j % p.num_pixels;
    double x2 = p.skewer_x[isk];
    double y2 = p.skewer_y[isk];
    double z2 = p.dz_pix * (iz + 0.5);

    // Separation vector components.
    double dx12 = x2 - x1;
    double dy12 = y2 - y1;
    double dz12 = z2 - z1;

    // The gaussian terms.
    double x_perp = sqrt(dx12*dx12 + dy12*dy12) / p.l_perp;
    double x_para = fabs(dz12) / p.l_para;

    // Look up gaussian product.
    double gp = p.gt->f(x_perp) * p.gt->f(x_para);
    // Don't forget the variance.
    double sn = p.sigma * gp;

    // Add noise if we're on the diag.
    if (i == j) {
        sn += p.noise[i];
    }

    return sn;
}
