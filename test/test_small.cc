/*
2014, the dachshund authors.

An idealized case to test reconstruction.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>

#include "dachshund.h"

int
main()
{
    //
    // Set params.
    //

    int num_skewers = 4;
    int num_pixels = 4;
    int nx = 16;
    int ny = 16;
    int nz = 16;
    double lx = 1.0;
    double ly = 1.0;
    double lz = 1.0;
    double sigma_F = 0.5;
    double l_perp = 0.125;
    double l_para = 0.25;
    double pcg_tol = 1.0e-3;
    int pcg_max_iter = 100;
    std::string output_path = "test_map.bin";

    int64_t num_cells = (int64_t)nx * (int64_t)ny * (int64_t)nz;
    int64_t num_pixel_elements = (int64_t)num_pixels * (int64_t)num_skewers;
    double dx = lx / nx;
    double dy = ly / ny;
    double dz = lz / nz;
    double dz_pix = lz / num_pixels;

    //
    // Set skewer data
    //

    double *skewer_x = new double[num_skewers];
    double *skewer_y = new double[num_skewers];
    double *data_vec = new double[num_pixel_elements];
    double *noise_vec = new double[num_pixel_elements];

    skewer_x[0] = 0.3;
    skewer_y[0] = 0.3;
    skewer_x[1] = 0.3;
    skewer_y[1] = 0.6;
    skewer_x[2] = 0.6;
    skewer_y[2] = 0.3;
    skewer_x[3] = 0.6;
    skewer_y[3] = 0.6;

    double dd[16] = {
        1.0, -1.0, 1.0, -1.0,
        -1.0, 1.0, -1.0, 1.0,
        1.0, -1.0, 1.0, -1.0,
        -1.0, 1.0, -1.0, 1.0
    };

    double nn[16] = {
        0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1
    };

    for (int i = 0; i < num_pixel_elements; ++i) {
        data_vec[i] = dd[i];
        noise_vec[i] = nn[i];
    }

    // Create gpt.
    int gpt_n = 100;
    double gpt_dx = 5.0 / gpt_n;
    GPT *gpt = new GPT(gpt_n, gpt_dx);

    //
    // Create matrix lookups.
    //

    double (*f)(int, int, void *);

    // S^MD matrix
    // params.
    SMDParams *p1 = new SMDParams;
    p1->nx = nx;
    p1->ny = ny;
    p1->nz = nz;
    p1->dx = dx;
    p1->dy = dy;
    p1->dz = dz;
    p1->ns = num_skewers;
    p1->np = num_pixels;
    p1->dz_pix = dz_pix;
    p1->skewer_x = skewer_x;
    p1->skewer_y = skewer_y;
    p1->l_perp = l_perp;
    p1->l_para = l_para;
    p1->sigma = sigma_F;
    p1->gpt = gpt;

    // Grab function pointer.
    f = &smd_element_func;

    LUM *smd = new LUM(num_cells, num_pixel_elements, f, p1);

    // S^DD + N matrix
    // params.
    SDDNParams *p2 = new SDDNParams;
    p2->ns = num_skewers;
    p2->np = num_pixels;
    p2->dz_pix = dz_pix;
    p2->skewer_x = skewer_x;
    p2->skewer_y = skewer_y;
    p2->l_perp = l_perp;
    p2->l_para = l_para;
    p2->sigma = sigma_F;
    p2->n = noise_vec;
    p2->gpt = gpt;

    // Grab function pointer.
    f = &sddn_element_func;

    LUM *sddn = new LUM(num_pixel_elements, num_pixel_elements, f, p2);

    // Allocate work arrays.
    double *a_j = new double[num_pixel_elements];
    double *d = new double[num_pixel_elements];
    double *r = new double[num_pixel_elements];
    double *s = new double[num_pixel_elements];
    double *q = new double[num_pixel_elements];
    // Allocate map.
    double *map_vec = new double[num_cells];

    // Zero out map.
    for (int i = 0; i < num_cells; ++i) {
        map_vec[i] = 0.0;
    }

    printf("Entering main loop with n_i = %i, n_j = %i.\n",
        (int)num_cells, (int)num_pixel_elements);

    // Set initial a_j guess.
    for (int i = 0; i < num_pixel_elements; ++i) {
        if (i == 0) {
            a_j[i] = 1.0 / (*sddn)(i, i);
        }
        else {
            a_j[i] = 0.0;
        }
    }

    // The main loop, first over cells, then pixels.
    for (int j = 0; j < num_pixel_elements; ++j) {

        // DEBUG
        printf("pixel %i -- %.2f %%\n", (int)j, 100.0*j/num_pixel_elements);

        // PCG for the j-th row of (S^DD + N)^-1.
        pcg_invert(j, sddn, a_j, pcg_max_iter, pcg_tol, d, r, s, q);

        // Dot with data vector.
        double b_j = 0.0;
        for (int k = 0; k < num_pixel_elements; ++k) {
            b_j += a_j[k] * data_vec[k];
        }

        // Add contribution to map.
        for (int i = 0; i < num_cells; ++i) {
            map_vec[i] += (*smd)(i, j) * b_j;
        }
    }

    // Write map field.
    printf("Writing map file %s.\n", output_path.c_str());
    FILE *outfile = fopen(output_path.c_str(), "wb");
    fwrite(map_vec, sizeof(double), num_cells, outfile);
    fclose(outfile);

    return 0;
}
