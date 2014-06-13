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
open_file_r(const std::string path, FILE *f)
{

    if (f == NULL) {
        printf("Could not open file %s\n", path.c_str());
        exit(1);
    }
}


void
read_data(double * const pixel_x, double * const pixel_y,
    double * const pixel_z,
    double * const pixel_data, double * const pixel_weights)
{
    FILE *infile = NULL;

    infile = fopen(p.pixel_x_path.c_str(), "r");
    if (infile == NULL) { printf("Could not open file %s\n", p.pixel_x_path.c_str()); exit(1); }
    fread(pixel_x, sizeof(double), p.num_pixels, infile);
    fclose(infile);

    infile = fopen(p.pixel_y_path.c_str(), "r");
    if (infile == NULL) { printf("Could not open file %s\n", p.pixel_y_path.c_str()); exit(1); }
    fread(pixel_y, sizeof(double), p.num_pixels, infile);
    fclose(infile);

    infile = fopen(p.pixel_z_path.c_str(), "r");
    if (infile == NULL) { printf("Could not open file %s\n", p.pixel_z_path.c_str()); exit(1); }
    fread(pixel_z, sizeof(double), p.num_pixels, infile);
    fclose(infile);

    infile = fopen(p.pixel_data_path.c_str(), "r");
    if (infile == NULL) { printf("Could not open file %s\n", p.pixel_data_path.c_str()); exit(1); }
    fread(pixel_data, sizeof(double), p.num_pixels, infile);
    fclose(infile);

    infile = fopen(p.pixel_weights_path.c_str(), "r");
    if (infile == NULL) { printf("Could not open file %s\n", p.pixel_weights_path.c_str()); exit(1); }
    fread(pixel_weights, sizeof(double), p.num_pixels, infile);
    fclose(infile);
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
    ds_params_print();

    //
    // Read in skewers.
    //

    double *pixel_x = new double[p.num_pixels];
    double *pixel_y = new double[p.num_pixels];
    double *pixel_z = new double[p.num_pixels];
    double *pixel_data = new double[p.num_pixels];
    double *pixel_weights = new double[p.num_pixels];

    puts("Reading skewer files.");
    read_data(pixel_x, pixel_y, pixel_z, pixel_data, pixel_weights);

    //
    // Loop setup.
    //

    // Apply weights to data...
    for (i = 0; i < p.num_pixels; ++i) {
        pixel_data[i] *= pixel_weights[i];
    }

    // Setup pixels.
    ds_lum_init(pixel_x, pixel_y, pixel_z, pixel_weights);

    delete [] pixel_x;
    delete [] pixel_y;
    delete [] pixel_z;
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
