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

void
read_data(double * const skewer_x, double * const skewer_y,
    double * const pixel_data, double * const pixel_weights)
{
    FILE *skewer_file = fopen(p.skewer_x_path.c_str(), "r");
    //check(skewer_file, "Could not open file %s.", p.skewer_x_path.c_str());
    fread(skewer_x, sizeof(double), p.num_skewers, skewer_file);
    fclose(skewer_file);

    skewer_file = fopen(p.skewer_y_path.c_str(), "r");
    //check(skewer_file, "Could not open file %s.", p.skewer_y_path.c_str());
    fread(skewer_y, sizeof(double), p.num_skewers, skewer_file);
    fclose(skewer_file);

    FILE *data_file = fopen(p.pixel_data_path.c_str(), "r");
    //check(data_file, "Could not open file %s.", p.pixel_data_path.c_str());
    fread(pixel_data, sizeof(double), p.pix_n, data_file);
    fclose(data_file);

    FILE *weights_file = fopen(p.pixel_weights_path.c_str(), "r");
    //check(weights_file, "Could not open file %s.", p.pixel_weights_path.c_str());
    fread(pixel_weights, sizeof(double), p.pix_n, weights_file);
    fclose(weights_file);
}

int
main(int argc, char **argv)
{
    int i;
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
    read_data(skewer_x, skewer_y, pixel_data, pixel_weights);

    //
    // Loop setup.
    //

    // Apply weights to data...
    for (i = 0; i < p.pix_n; ++i) {
        pixel_data[i] *= pixel_weights[i];
    }

    // Setup pixels.
    ds_lum_init(skewer_x, skewer_y, pixel_weights);

    delete [] skewer_x;
    delete [] skewer_y;
    delete [] pixel_weights;

    // Send off to compute the map.
    double *map = new double[p.map_n];
    ds_compute_map(pixel_data, map);

    // Write map field.
    printf("Writing map file %s.\n", p.map_path.c_str());
    FILE *map_file = fopen(p.map_path.c_str(), "w");
    fwrite(map, sizeof(double), p.map_n, map_file);
    fclose(map_file);

    if (p.option_noise_covar) {
        puts("Computing covariance in N >> S limit.");
        ds_smp_w_spm(map);
        map_file = fopen("map_covar_nd.bin", "w");
        fwrite(map, sizeof(double), p.map_n, map_file);
        fclose(map_file);
    }

    if (p.option_map_covar) {
        puts("Computing covariance diag.");
        ds_map_covar_diag(map);
        map_file = fopen("map_covar.bin", "w");
        fwrite(map, sizeof(double), p.map_n, map_file);
        fclose(map_file);
    }

    printf("[TIME] %g ms total.\n", total_timer->elapsed());

    return 0;
}
