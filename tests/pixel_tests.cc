/*
2014, the dachshund authors.
*/

#include <cmath>

#include "catch.hpp"

#include "dachshund.h"
#include "test_utils.h"

const int n = 10;

TEST_CASE("Pixel conversions", "[pixel]") {
    double *nvec = new double[n];
    double *dvec = new double[n];
    Pixel *pixels = new Pixel[n];

    for (int i = 0; i < n; ++i) {
        nvec[i] = rng();
        dvec[i] = rng();

        pixels[i].x = 0.0;
        pixels[i].y = 0.0;
        pixels[i].z = 0.0;
        pixels[i].n = nvec[i];
        pixels[i].d = dvec[i];
    }

    SECTION("pixel to npixel") {
        NPixel *npixels = (NPixel *) pixels;
        pixel_to_npixel(n, pixels, npixels);

        for (int i = 0; i < n; ++i) {
            REQUIRE( npixels[i].N == nvec[i] * nvec[i] );
            REQUIRE( npixels[i].d == dvec[i] );
        }

        npixel_to_pixel(n, npixels, pixels);

        for (int i = 0; i < n; ++i) {
            REQUIRE( pixels[i].n == nvec[i] );
            REQUIRE( pixels[i].d == dvec[i] );
        }
    }

    SECTION("pixel to wpixel") {
        WPixel *wpixels = (WPixel *) pixels;
        pixel_to_wpixel(n, pixels, wpixels);

        for (int i = 0; i < n; ++i) {
            REQUIRE( wpixels[i].w == 1.0 / nvec[i] );
            REQUIRE( wpixels[i].wd == dvec[i] / nvec[i] );
        }

        wpixel_to_pixel(n, wpixels, pixels);

        for (int i = 0; i < n; ++i) {
            // add some tol.
            REQUIRE( dtol(pixels[i].n, nvec[i], 0.0, 1e-15) );
            REQUIRE( dtol(pixels[i].d, dvec[i], 0.0, 1e-15) );
        }
    }

    delete [] nvec;
    delete [] dvec;
    delete [] pixels;
}
