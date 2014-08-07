/*
2014, the dachshund authors.
*/

#ifndef __DS_PIXEL_H__
#define __DS_PIXEL_H__

struct Point {
    double x;
    double y;
    double z;
};

// the default pixel

struct Pixel {
    double x;
    double y;
    double z;
    double n;
    double d;
};

Point *
points_from_pixels_alloc(const int n, Pixel * const pixels);

// alternate pixel formats for performance in loops

struct NPixel {
    double x;
    double y;
    double z;
    double N;
    double d;
};

void
pixel_to_npixel(const int n, Pixel * const p, NPixel * const np);

void
npixel_to_pixel(const int n, NPixel * const np, Pixel * const p);

struct WPixel{
    double x;
    double y;
    double z;
    double w;
    double wd;
};

void
pixel_to_wpixel(const int n, Pixel * const p, WPixel * const wp);

void
wpixel_to_pixel(const int n, WPixel * const wp, Pixel * const p);

#endif
