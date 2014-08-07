/*
2014, the dachshund authors.
*/

#include <cmath>

#include "pixel.h"

Point *
points_from_pixels_alloc(const int n, Pixel * const pixels)
{
    Point * const points = new Point[n];
    for (int i = 0; i < n; ++i) {
        points[i].x = pixels[i].x;
        points[i].y = pixels[i].y;
        points[i].z = pixels[i].z;
    }
    return points;
}

// pixel formats...

void
pixel_to_npixel(const int n, Pixel * const p, NPixel * const np)
{
    for (int i = 0; i < n; ++i) {
        np[i].x = p[i].x;
        np[i].y = p[i].y;
        np[i].z = p[i].z;
        np[i].N = p[i].n * p[i].n;
        np[i].d = p[i].d;
    }
}

void
npixel_to_pixel(const int n, NPixel * const np, Pixel * const p)
{
    for (int i = 0; i < n; ++i) {
        p[i].x = np[i].x;
        p[i].y = np[i].y;
        p[i].z = np[i].z;
        p[i].n = sqrt(np[i].N);
        p[i].d = np[i].d;
    }
}

void
pixel_to_wpixel(const int n, Pixel * const p, WPixel * const wp)
{
    for (int i = 0; i < n; ++i) {
        wp[i].x = p[i].x;
        wp[i].y = p[i].y;
        wp[i].z = p[i].z;
        double ni = p[i].n;
        wp[i].w = 1.0 / ni;
        wp[i].wd = p[i].d / ni;
    }
}

void
wpixel_to_pixel(const int n, WPixel * const wp, Pixel * const p)
{
    for (int i = 0; i < n; ++i) {
        p[i].x = wp[i].x;
        p[i].y = wp[i].y;
        p[i].z = wp[i].z;
        double wi = wp[i].w;
        p[i].n = 1.0 / wi;
        p[i].d = wp[i].wd / wi;
    }
}

