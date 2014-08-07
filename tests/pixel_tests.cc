/*
2014, the dachshund authors.
*/

#include <cmath>

#include "dachshund.h"
#include "test_utils.h"

const int n = 10;

struct RandomPixels
{
    double *nvec;
    double *dvec;
    Pixel *pixels;

    RandomPixels() {
        nvec = new double[n];
        dvec = new double[n];
        pixels = new Pixel[n];

        for (int i = 0; i < n; ++i) {
            nvec[i] = rng();
            dvec[i] = rng();

            pixels[i].x = 0.0;
            pixels[i].y = 0.0;
            pixels[i].z = 0.0;
            pixels[i].n = nvec[i];
            pixels[i].d = dvec[i];
        }
    }
    ~RandomPixels() {
        delete [] nvec;
        delete [] dvec;
        delete [] pixels;
    }
};

TEST_FIXTURE(RandomPixels, pixel_npixel_same)
{
    NPixel *npixels = (NPixel *) pixels;

    pixel_to_npixel(n, pixels, npixels);

    for (int i = 0; i < n; ++i) {
        CHECK_EQUAL(npixels[i].N, nvec[i] * nvec[i]);
        CHECK_EQUAL(npixels[i].d, dvec[i]);
    }

    npixel_to_pixel(n, npixels, pixels);

    for (int i = 0; i < n; ++i) {
        CHECK_EQUAL(pixels[i].n, nvec[i]);
        CHECK_EQUAL(pixels[i].d, dvec[i]);
    }
}

TEST_FIXTURE(RandomPixels, pixel_npixel_copy)
{
    NPixel *npixels = new NPixel[n];

    pixel_to_npixel(n, pixels, npixels);

    for (int i = 0; i < n; ++i) {
        CHECK_EQUAL(npixels[i].N, nvec[i] * nvec[i]);
        CHECK_EQUAL(npixels[i].d, dvec[i]);
    }

    npixel_to_pixel(n, npixels, pixels);

    for (int i = 0; i < n; ++i) {
        CHECK_EQUAL(pixels[i].n, nvec[i]);
        CHECK_EQUAL(pixels[i].d, dvec[i]);
    }
}

TEST_FIXTURE(RandomPixels, pixel_wpixel_same)
{
    WPixel *wpixels = (WPixel *) pixels;

    pixel_to_wpixel(n, pixels, wpixels);

    for (int i = 0; i < n; ++i) {
        CHECK_EQUAL(wpixels[i].w, 1.0 / nvec[i]);
        CHECK_EQUAL(wpixels[i].wd, dvec[i] / nvec[i]);
    }

    wpixel_to_pixel(n, wpixels, pixels);

    for (int i = 0; i < n; ++i) {
        CHECK_DTOL(pixels[i].n, nvec[i], 0.0, 1.0e-15);
        CHECK_DTOL(pixels[i].d, dvec[i], 0.0, 1.0e-15);
    }
}

TEST_FIXTURE(RandomPixels, pixel_wpixel_copy)
{
    WPixel *wpixels = (WPixel *) pixels;

    pixel_to_wpixel(n, pixels, wpixels);

    for (int i = 0; i < n; ++i) {
        CHECK_EQUAL(wpixels[i].w, 1.0 / nvec[i]);
        CHECK_EQUAL(wpixels[i].wd, dvec[i] / nvec[i]);
    }

    wpixel_to_pixel(n, wpixels, pixels);

    for (int i = 0; i < n; ++i) {
        CHECK_DTOL(pixels[i].n, nvec[i], 0.0, 1.0e-15);
        CHECK_DTOL(pixels[i].d, dvec[i], 0.0, 1.0e-15);
    }
}
