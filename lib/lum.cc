/*
2014, the dachshund authors.
*/

#include <cstdio>
#include <cstdlib>
#include <math.h>

#include "params.h"

#include "lum.h"

// global static variables
static DSFuncTable gt;
static DSPixel *pixels;
static DSPoint *map_coords;

void
ds_lum_gt_init(const int n, const double dx)
{
    gt.n = n;
    gt.dx = dx;
    gt.table = new double[n];
    for (int i = 0; i < n; ++i) {
        double x = i * dx;
        gt.table[i] = exp(-x);
    }
}

void
ds_lum_gt_free()
{
    delete gt.table;
}

void
ds_lum_pixels_init(const int num_skewers,
    const double * const skewer_x, const double * const skewer_y,
    const int num_pixels, const double pix_dz,
    const double * const pixel_weights)
{
    int pix_n = num_skewers * num_pixels;
    pixels = new DSPixel[pix_n];

    for (int i = 0; i < pix_n; ++i) {
        int isk = i / num_pixels;
        int iz = i % num_pixels;
        pixels[i].x = skewer_x[isk];
        pixels[i].y = skewer_y[isk];
        pixels[i].z = pix_dz * (iz + 0.5);
        pixels[i].w = pixel_weights[i];
    }
}

void
ds_lum_pixels_free()
{
    delete pixels;
}

void
ds_lum_map_coords_init(const int map_nx, const int map_ny, const int map_nz,
    const int map_dx, const int map_dy, const int map_dz)
{
    int map_n = map_nx * map_ny * map_nz;
    map_coords = new DSPoint[map_n];

    for (int ix = 0; ix < map_nx; ++ix) {
        for (int iy = 0; iy < map_ny; ++iy) {
            for (int iz = 0; iz < map_nz; ++iz) {
                int i = (ix * map_ny + iy) * map_nz + iz;
                map_coords[i].x = map_dx * (ix + 0.5);
                map_coords[i].y = map_dy * (iy + 0.5);
                map_coords[i].z = map_dz * (iz + 0.5);
            }
        }
    }
}

void
ds_lum_map_coords_free()
{
    delete map_coords;
}

/*
Linear interpolation, given x and table data.

We assume x >= 0. If x is larger than the last tabulated point, we return 0.

Parameters
x : the eval point
nx : the number of samples
dx: the sample point separation
f_table : array containing the function value at sample points.
*/
static inline double
ds_linterp(const double x, const int nx, const double dx,
    const double * const f_table)
{
    // in index space.
    double xx = x / dx;
    // decimate to get low edge.
    int ix = xx;
    // Check if we are outside table range.
    if (ix > nx - 2) {
        return 0.0;
    }
    // the linterp expression.
    return f_table[ix] * (ix + 1 - xx) + f_table[ix+1] * (xx - ix);
}


double
wsppi_lookup(const int i, const int j)
{
    static double dx, dy, dz;
    static double x_perp_2, x_para_2;
    static double g_perp, g_para, sn;

    // Separation vector components.
    dx = pixels[i].x - pixels[j].x;
    dy = pixels[i].y - pixels[j].y;
    dz = pixels[i].z - pixels[j].z;

    // The gaussian terms.
    x_perp_2 = (dx*dx + dy*dy) / p.corr_l_perp_2;
    x_para_2 = dz*dz / p.corr_l_para_2;

    // S = sigma^2 exp(...) exp(...)
    g_perp = ds_linterp(x_perp_2, gt.n, gt.dx, gt.table);
    g_para = ds_linterp(x_para_2, gt.n, gt.dx, gt.table);
    sn = p.corr_var_s * g_perp * g_para;

    // Add weight and identity.
    if (i == j) {
        sn = pixels[i].w * sn + 1.0;
    }

    return sn;
}

double
smp_lookup(const int i, const int j)
{
    static double dx, dy, dz;
    static double x_perp_2, x_para_2;
    static double g_perp, g_para;

    // Separation vector components.
    dx = map_coords[i].x - pixels[j].x;
    dy = map_coords[i].y - pixels[j].y;
    dz = map_coords[i].z - pixels[j].z;

    // The gaussian terms.
    x_perp_2 = dx*dx + dy*dy / p.corr_l_perp_2;
    x_para_2 = dz*dz / p.corr_l_para_2;

    // S = sigma^2 exp(...) exp(...)
    g_perp = ds_linterp(x_perp_2, gt.n, gt.dx, gt.table);
    g_para = ds_linterp(x_para_2, gt.n, gt.dx, gt.table);

    return p.corr_var_s * g_perp * g_para;
}

double
wspm_lookup(const int i, const int j)
{
    static double dx, dy, dz;
    static double x_perp_2, x_para_2;
    static double g_perp, g_para;

    // Separation vector components.
    dx = pixels[i].x - map_coords[j].x;
    dy = pixels[i].y - map_coords[j].y;
    dz = pixels[i].z - map_coords[j].z;

    // The gaussian terms.
    x_perp_2 = dx*dx + dy*dy / p.corr_l_perp_2;
    x_para_2 = dz*dz / p.corr_l_para_2;

    // S = sigma^2 exp(...) exp(...)
    g_perp = ds_linterp(x_perp_2, gt.n, gt.dx, gt.table);
    g_para = ds_linterp(x_para_2, gt.n, gt.dx, gt.table);

    return pixels[i].w * p.corr_var_s * g_perp * g_para;
}
