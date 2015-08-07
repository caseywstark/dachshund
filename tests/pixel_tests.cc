/*
2014, the dachshund authors.
*/

#include "gtest/gtest.h"

#include "dachshund.h"
#include "test_utils.h"

class PixelTest : public testing::Test {
  protected:
    const static int n = 10;

    virtual void SetUp() {
        n_vec.resize(n);
        d_vec.resize(n);
        pixels = new Pixel[n];

        for (int i = 0; i < n; ++i) {
            n_vec[i] = rng();
            d_vec[i] = rng();

            pixels[i].x = 0.0;
            pixels[i].y = 0.0;
            pixels[i].z = 0.0;
            pixels[i].n = n_vec[i];
            pixels[i].d = d_vec[i];
        }
    }

    virtual void TearDown() {
        delete [] pixels;
    }

    std::vector<double> n_vec;
    std::vector<double> d_vec;
    Pixel *pixels;
};

TEST_F(PixelTest, NPixel_conversion) {
    // Reuse the same block.
    NPixel *npixels = (NPixel *) pixels;
    pixel_to_npixel(n, pixels, npixels);

    for (int i = 0; i < n; ++i) {
        EXPECT_EQ(n_vec[i] * n_vec[i], npixels[i].N);
        EXPECT_EQ(d_vec[i], npixels[i].d);
    }

    npixel_to_pixel(n, npixels, pixels);

    for (int i = 0; i < n; ++i) {
        EXPECT_EQ(n_vec[i], pixels[i].n);
        EXPECT_EQ(d_vec[i], pixels[i].d);
    }
}

TEST_F(PixelTest, WPixel_conversion) {
    WPixel *wpixels = (WPixel *) pixels;
    pixel_to_wpixel(n, pixels, wpixels);

    for (int i = 0; i < n; ++i) {
        EXPECT_EQ(1.0 / n_vec[i], wpixels[i].w);
        EXPECT_EQ(d_vec[i] / n_vec[i], wpixels[i].wd);
    }

    wpixel_to_pixel(n, wpixels, pixels);

    for (int i = 0; i < n; ++i) {
        // Add some tol.
        EXPECT_EQTOL(n_vec[i], pixels[i].n, 0.0, 1e-15);
        EXPECT_EQTOL(d_vec[i], pixels[i].d, 0.0, 1e-15);
    }
}
