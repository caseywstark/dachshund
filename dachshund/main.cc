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
    double *pixel_data = new double[p.pix_n];
    double *pixel_weights = new double[p.pix_n];

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
    fread(pixel_data, sizeof(double), p.pix_n, data_file);
    fclose(data_file);

    // Skewer weights
    FILE *weights_file = fopen(p.pixel_weights_path.c_str(), "r");
    if (weights_file == NULL) {
        fprintf(stderr, "Could not load file %s.\n", p.pixel_weights_path.c_str());
        exit(1);
    }
    fread(pixel_weights, sizeof(double), p.pix_n, weights_file);
    fclose(weights_file);

    //
    // Loop setup.
    //

    // Apply weights to data...
    for (i = 0; i < p.pix_n; ++i) {
        pixel_data[i] *= pixel_weights[i];
    }

    // Setup pixels.
    ds_lum_pixels_init(p.num_skewers, skewer_x, skewer_y, p.num_pixels,
        p.pix_dz, pixel_weights);
    ds_lum_map_coords_init(p.map_nx, p.map_ny, p.map_nz, p.map_dx, p.map_dy,
        p.map_dz);

    delete [] skewer_x;
    delete [] skewer_y;
    delete [] pixel_weights;

    int gt_n = 50;
    double gt_dx = 2.0 / gt_n;
    ds_lum_gt_init(gt_n, gt_dx);

    //
    // PCG step
    //

    printf("PCG with pix_n = %i.\n", p.pix_n);

    double *x = new double[p.pix_n];
    double *b = pixel_data;

    for (i = 0; i < p.pix_n; ++i) {
        x[i] = 0.0;
    }

    // We care about timing starting here.
    Timer *pcg_timer = new Timer();

    pcg(p.pix_n, &wsppi_lookup, x, b, p.pcg_max_iter, p.pcg_tol, true);

    printf("PCG time: %g ms.\n", pcg_timer->elapsed());

    //
    // Matrix multiply for map.
    //

    // Allocate map.
    double *map = new double[p.map_n];

    printf("Multiply for map_n = %i, total n = %.2e.\n", p.map_n, (double)p.map_n * p.pix_n);

    // Each thread gets a block of cells.
#if defined(_OPENMP)
    #pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < p.map_n; ++i) {
        map[i] = 0.0;
        for (j = 0; j < p.pix_n; ++j) {
            map[i] += smp_lookup(i, j) * x[j];
        }
    }

    // Write map field.
    printf("Writing map file %s.\n", p.map_path.c_str());
    FILE *map_file = fopen(p.map_path.c_str(), "w");
    fwrite(map, sizeof(double), p.map_n, map_file);
    fclose(map_file);

    if (p.option_compute_covar) {
        puts("Computing covariance diag.");
        for (int i = 0; i < p.map_n; ++i) {
            printf("cell %i.", i);
            // Compute x_i = (S +N)^{-1}_{ij} S^{pm}_{j beta}
            pcg_covar(p.pix_n, &wsppi_lookup, x, &wspm_lookup, i,
                p.pcg_max_iter, p.pcg_tol, false);

            // Next step C_{alpha i}  = S^{mp}_{alpha i} x_i
            // Store in the map vector.
            map[i] = 0.0;
#if defined(_OPENMP)
            #pragma omp parallel for private(j)
#endif
            for (int j = 0; j < p.pix_n; ++j) {
                map[i] += smp_lookup(i, j) * x[j];
            }
        }

        map_file = fopen("map_covar.bin", "w");
        fwrite(map, sizeof(double), p.map_n, map_file);
        fclose(map_file);
    }

    printf("Total time: %g ms.\n", total_timer->elapsed());

    return 0;
}
