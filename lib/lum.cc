/*
2014, the dachshund authors.
*/

#include <cstdio>
#include <math.h>

#include "params.h"

#include "lum.h"

double
smp_lookup(const int i, const int j)
{
    // Separation vector components.
    double dx = p.map[i].x - p.pixels[j].x;
    double dy = p.map[i].y - p.pixels[j].y;
    double dz = p.map[i].z - p.pixels[j].z;

    // The gaussian terms.
    double x_perp = sqrt(dx*dx + dy*dy) / p.l_perp;
    double x_para = fabs(dz) / p.l_para;

    // S = sigma^2 exp(...) exp(...)
    double sn = p.sigma * p.gt->f(x_perp) * p.gt->f(x_para);
    return sn;
}

double
wsppi_lookup(const int i, const int j)
{
    // Separation vector components.
    double dx = p.pixels[i].x - p.pixels[j].x;
    double dy = p.pixels[i].y - p.pixels[j].y;
    double dz = p.pixels[i].z - p.pixels[j].z;

    // The gaussian terms.
    double x_perp = sqrt(dx*dx + dy*dy) / p.l_perp;
    double x_para = fabs(dz) / p.l_para;

    // S = sigma^2 exp(...) exp(...)
    double sn = p.sigma * p.gt->f(x_perp) * p.gt->f(x_para);

    // Add weight and identity.
    if (i == j) {
        sn = p.pixels[i].w * sn + 1.0;
    }

    return sn;
}

double
wspm_lookup(const int i, const int j)
{
    // Separation vector components.
    double dx = p.pixels[i].x - p.map[j].x;
    double dy = p.pixels[i].y - p.map[j].y;
    double dz = p.pixels[i].z - p.map[j].z;

    // The gaussian terms.
    double x_perp = sqrt(dx*dx + dy*dy) / p.l_perp;
    double x_para = fabs(dz) / p.l_para;

    // S = sigma^2 exp(...) exp(...)
    double sn = p.pixels[i].w * p.sigma * p.gt->f(x_perp) * p.gt->f(x_para);
    return sn;
}
