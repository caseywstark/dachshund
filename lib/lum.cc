/*
2014, the dachshund authors.
*/

#include <cstdio>
#include <math.h>

#include "lum.h"

LUM::LUM(int n, int m, double (*e)(int i, int j, void *params), void *p)
{
    this->n = n;
    this->m = m;
    element_func = e;
    params = p;
}

double
LUM::operator()(int i, int j) const
{
    double e = (*element_func)(i, j, params);
    return e;
}

double
smd_element_func(int i, int j, void *params)
{
    // unpack params.
    SMDParams *p = (SMDParams *) params;
    int nx = p->nx;
    int ny = p->ny;
    int nz = p->nz;
    double dx = p->dx;
    double dy = p->dy;
    double dz = p->dz;
    int ns = p->ns;
    int np = p->np;
    double dz_pix = p->dz_pix;
    double *skewer_x = p->skewer_x;
    double *skewer_y = p->skewer_y;
    double l_perp = p->l_perp;
    double l_para = p->l_para;
    double sigma = p->sigma;
    GT *gt = p->gt;

    // grid flat index i -> grid ix, iy, iz
    // probably a faster way to do this...
    int iz = i % nz;
    int iy = (i - iz) / nz % ny;
    int ix = (i - iz - iy * nz) / (ny * nz) % nx;

    double x1 = dx * (ix + 0.5);
    double y1 = dy * (iy + 0.5);
    double z1 = dz * (iz + 0.5);

    int isk = j / np;
    int iz2 = j % np;
    double x2 = skewer_x[isk];
    double y2 = skewer_y[isk];
    double z2 = dz_pix * (iz2 + 0.5);

    // Separation vector components.
    double dx12 = x2 - x1;
    double dy12 = y2 - y1;
    double dz12 = z2 - z1;

    // The gaussian terms.
    double x_perp = sqrt(dx12*dx12 + dy12*dy12) / l_perp;
    double x_para = fabs(dz12) / l_para;

    // Look up gaussian product.
    double gp = (*gt)(x_perp) * (*gt)(x_para);
    // Don't forget the variance.
    double sn = sigma * gp;

    return sn;
}


double
sddn_element_func(int i, int j, void *params)
{
    // unpack params.
    SDDNParams *p = (SDDNParams *) params;
    int ns = p->ns;
    int np = p->np;
    double dz_pix = p->dz_pix;
    double *skewer_x = p->skewer_x;
    double *skewer_y = p->skewer_y;
    double l_perp = p->l_perp;
    double l_para = p->l_para;
    double sigma = p->sigma;
    double *n = p->n;
    GT *gt = p->gt;

    int isk, iz;
    isk = i / np;
    iz = i % np;
    double x1 = skewer_x[isk];
    double y1 = skewer_y[isk];
    double z1 = dz_pix * (iz + 0.5);

    isk = j / np;
    iz = j % np;
    double x2 = skewer_x[isk];
    double y2 = skewer_y[isk];
    double z2 = dz_pix * (iz + 0.5);

    // Separation vector components.
    double dx12 = x2 - x1;
    double dy12 = y2 - y1;
    double dz12 = z2 - z1;

    // The gaussian terms.
    double x_perp = sqrt(dx12*dx12 + dy12*dy12) / l_perp;
    double x_para = fabs(dz12) / l_para;

    // Look up gaussian product.
    double gp = (*gt)(x_perp) * (*gt)(x_para);
    // Don't forget the variance.
    double sn = sigma * gp;

    // Add noise if we're on the diag.
    if (i == j) {
        sn += n[i];
    }

    return sn;
}
