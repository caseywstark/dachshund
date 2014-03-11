/*
2014, the dachshund authors.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>

#include <omp.h>

#include "dachshund.h"

void
print_usage()
{
    puts("Usage: dachshund.ex input");
    puts("    input     path to config file.");
}

int
main(int argc, char **argv)
{
    int i, j;
    Timer *total_timer = new Timer();

    // Fancy command line parsing.
    if (argc != 2) {
        fputs("Bad number of arguments. Expected 1.\n", stderr);
        print_usage();
        exit(1);
    }

    std::string config_path = argv[1];

    //
    // Read config file and set params.
    //

    ds_params_init(config_path);

    //
    // Read in skewers.
    //

    double *skewer_x = new double[p.num_skewers];
    double *skewer_y = new double[p.num_skewers];
    double *pixel_data = new double[p.num_pixel_elements];
    double *pixel_weights = new double[p.num_pixel_elements];

    puts("Reading skewer files.");

    // Skewer positions
    FILE *skewer_file = fopen(p.skewer_x_path.c_str(), "r");
    if (skewer_file == NULL) {
        fprintf(stderr, "Could not load file %s.\n", p.skewer_x_path.c_str());
        exit(1);
    }
    fread(skewer_x, sizeof(double), p.num_skewers, skewer_file);
    fclose(skewer_file);

    skewer_file = fopen(p.skewer_y_path.c_str(), "r");
    if (skewer_file == NULL) {
        fprintf(stderr, "Could not load file %s.\n", p.skewer_y_path.c_str());
        exit(1);
    }
    fread(skewer_y, sizeof(double), p.num_skewers, skewer_file);
    fclose(skewer_file);

    // Skewer fluxes
    FILE *data_file = fopen(p.pixel_data_path.c_str(), "r");
    if (data_file == NULL) {
        fprintf(stderr, "Could not load file %s.\n", p.pixel_data_path.c_str());
        exit(1);
    }
    fread(pixel_data, sizeof(double), p.num_pixel_elements, data_file);
    fclose(data_file);

    // Skewer weights
    FILE *weights_file = fopen(p.pixel_weights_path.c_str(), "r");
    if (weights_file == NULL) {
        fprintf(stderr, "Could not load file %s.\n", p.pixel_weights_path.c_str());
        exit(1);
    }
    fread(pixel_weights, sizeof(double), p.num_pixel_elements, weights_file);
    fclose(weights_file);

    //
    // Loop setup.
    //

    // Apply weights to data...
    for (int64_t i = 0; i < p.num_pixel_elements; ++i) {
        pixel_data[i] *= pixel_weights[i];
    }

    // Setup pixels.
    p.pixels = new DSPixel[p.num_pixel_elements];

    for (i = 0; i < p.num_pixel_elements; ++i) {
        int isk = i / p.num_pixels;
        int iz = i % p.num_pixels;
        p.pixels[i].x = skewer_x[isk];
        p.pixels[i].y = skewer_y[isk];
        p.pixels[i].z = p.dz_pix * (iz + 0.5);
        p.pixels[i].w = pixel_weights[i];
    }

    delete [] skewer_x;
    delete [] skewer_y;
    delete [] pixel_weights;

    // Setup cells.
    p.map = new DSPoint[p.num_cells];

    for (int ix = 0; ix < p.nx; ++ix) {
        for (int iy = 0; iy < p.ny; ++iy) {
            for (int iz = 0; iz < p.nz; ++iz) {
                int i = (ix * p.ny + iy) * p.nz + iz;
                p.map[i].x = p.dx * (ix + 0.5);
                p.map[i].y = p.dy * (iy + 0.5);
                p.map[i].z = p.dz * (iz + 0.5);
            }
        }
    }

    // Create gpt.
    int gt_n = 50;
    double gt_dx = 5.0 / gt_n;
    GT *gt = new GT(gt_n, gt_dx);
    p.gt = gt;

    // Correct PCG tolerance for the current number of pixels...
    double res_pcg_tol = p.pcg_tol / p.num_pixel_elements;

    printf("Entering inversion loop, num_pixel_elements = %i.\n",
        p.num_pixel_elements);

    double *x = new double[p.num_pixel_elements];
    double *b = pixel_data;

    for (i = 0; i < p.num_pixel_elements; ++i) {
        x[i] = 0.0;
    }

    // We care about timing starting here.
    Timer *loop_timer = new Timer();

    pcg(p.num_pixel_elements, &wsppi_lookup, x, b, p.pcg_max_iter, res_pcg_tol);

    printf("Main loop time: %g ms.\n", loop_timer->elapsed());

    //
    // Second loop to compute S b
    //

    // Allocate map.
    double *map = new double[p.num_cells];

    // Each thread gets a block of cells.
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < p.num_cells; ++i) {
        map[i] = 0.0;
        for (j = 0; j < p.num_pixel_elements; ++j) {
            map[i] += smp_lookup(i, j) * x[j];
        }
    }

    // Write map field.
    printf("Writing map file %s.\n", p.map_path.c_str());
    FILE *map_file = fopen(p.map_path.c_str(), "w");
    fwrite(map, sizeof(double), p.num_cells, map_file);
    fclose(map_file);

    if (p.compute_covar) {
        puts("Computing covariance diag.");
        for (int i = 0; i < p.num_cells; ++i) {
            // Compute x_i = (S +N)^{-1}_{ij} S^{pm}_{j beta}
            pcg_covar(p.num_pixel_elements, &wsppi_lookup, x, &wspm_lookup, i,
                p.pcg_max_iter, res_pcg_tol);

            // Next step C_{alpha i}  = S^{mp}_{alpha i} x_i
            // Store in the map vector.
            map[i] = 0.0;
#if defined(_OPENMP)
            #pragma omp parallel for private(j)
#endif
            for (int j = 0; j < p.num_pixel_elements; ++j) {
                map[i] += smp_lookup(i, j) * x[j];
            }
        }

        map_file = fopen("map_covar.bin", "w");
        fwrite(map, sizeof(double), p.num_cells, map_file);
        fclose(map_file);
    }

    printf("Total time: %g ms.\n", total_timer->elapsed());

    return 0;
}
