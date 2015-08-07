/*
2014, the dachshund authors.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>

#include <iostream>

#include "dachshund.h"

//
// Parameters, with default values.
//

// domain
double lx = 10.0;
double ly = 10.0;
double lz = 10.0;

int num_pixels = 100;

// map
int map_nx = 10;
int map_ny = 10;
int map_nz = 10;
int num_map_points = 1000;

// correlation
double corr_var_s = 1.0;
double corr_l_perp = 1.0;
double corr_l_para = 1.0;

// WF x solve
bool solve_via_pcg = true;

// PCG
int pcg_max_iter = 100;
double pcg_tol = 1.0;
int pcg_step_r = 50;

// options
bool option_map_covar = 0;
bool option_noise_covar = 0;

// paths
std::string pixel_data_path = "pixel_data.bin";
std::string map_path = "map.bin";

// end parameters

void
print_usage()
{
    puts("Usage: dachshund.ex input");
    puts("    input     path to config file.");
}

void
params_init(const std::string config_path)
{
    //
    // Parse input/config.
    //

    printf("Reading config file %s.\n", config_path.c_str());

    const static int max_line_length = 1000;
    char config_buf[max_line_length];
    char key[100];
    char value[100];

    // Open input file.
    FILE *config_file = fopen(config_path.c_str(), "r");
    // Make sure the input file is ok.
    if (config_file == NULL) {
        fprintf(stderr, "Could not load file %s.\n", config_path.c_str());
        exit(1);
    }

    // Have we reached the end of the file?
    while (!feof(config_file)) {
        // read current line.
        fgets(config_buf, max_line_length, config_file);

        // Handle blank or comment line.
        if (config_buf[0] == '\0' || config_buf[0] == '#') {
            continue;
        }

        // Parse normal line.
        sscanf(config_buf, "%s = %s", key, value);

        if (strcmp(key, "lx") == 0) {
            lx = atof(value);
        }
        else if (strcmp(key, "ly") == 0) {
            ly = atof(value);
        }
        else if (strcmp(key, "lz") == 0) {
            lz = atof(value);
        }
        else if (strcmp(key, "num_pixels") == 0) {
            num_pixels = atoi(value);
        }
        else if (strcmp(key, "map_nx") == 0) {
            map_nx = atoi(value);
        }
        else if (strcmp(key, "map_ny") == 0) {
            map_ny = atoi(value);
        }
        else if (strcmp(key, "map_nz") == 0) {
            map_nz = atoi(value);
        }
        else if (strcmp(key, "corr_var_s") == 0) {
            corr_var_s = atof(value);
        }
        else if (strcmp(key, "corr_l_para") == 0) {
            corr_l_para = atof(value);
        }
        else if (strcmp(key, "corr_l_perp") == 0) {
            corr_l_perp = atof(value);
        }
        else if (strcmp(key, "solve_via_pcg") == 0) {
            solve_via_pcg = atoi(value);
        }
        else if (strcmp(key, "pcg_tol") == 0) {
            pcg_tol = atof(value);
        }
        else if (strcmp(key, "pcg_max_iter") == 0) {
            pcg_max_iter = atoi(value);
        }
        else if (strcmp(key, "pcg_step_r") == 0) {
            pcg_step_r = atoi(value);
        }
        else if (strcmp(key, "option_noise_covar") == 0) {
            option_noise_covar = atoi(value);
        }
        else if (strcmp(key, "option_map_covar") == 0) {
            option_map_covar = atoi(value);
        }
        else if (strcmp(key, "pixel_data_path") == 0) {
            pixel_data_path = value;
        }
        else if (strcmp(key, "map_path") == 0) {
            map_path = value;
        }
        else {
            fprintf(stderr, "Unknown key '%s' in config.\n", key);
            exit(1);
        }
    }

    fclose(config_file);

    num_map_points = map_nx * map_ny * map_nz;
}

void
params_print()
{
    printf("\n");
    printf("Running dachshund with:\n");
    printf("  pixels/map:\n");
    printf("    lx ly lz = %f %f %f\n", lx, ly, lz);
    printf("    npix = %i\n", num_pixels);
    printf("    map nx ny nz = %i %i %i, n = %i\n", map_nx, map_ny, map_nz, num_map_points);
    printf("  signal:\n");
    printf("    corr var_s = %f\n", corr_var_s);
    printf("    l_perp = %f    l_para = %f\n", corr_l_perp, corr_l_para);
    printf("  WF x solve: %i\n", solve_via_pcg);
    printf("  PCG:\n");
    printf("    max iter = %i    step r = %i\n", pcg_max_iter, pcg_step_r);
    printf("    tol = %f\n", pcg_tol);
    printf("  options:\n");
    printf("    map covar = %i\n", option_map_covar);
    printf("  paths:\n");
    printf("    pixel data = %s\n", pixel_data_path.c_str());
    printf("    map = %s\n", map_path.c_str());
    printf("\n");
}


int
main(int argc, char **argv)
{
    //
    // Fancy command line parsing.
    //

    if (argc != 2) {
        fputs("Bad number of arguments. Expected 1.\n", stderr);
        print_usage();
        exit(1);
    }

    std::string config_path = argv[1];

    params_init(config_path);
    params_print();

    //
    // Reconstruction setup.
    //

    // Allocate everything we will pass to reconstruction here.
    Pixel * const pixels = new Pixel[num_pixels];
    Point * const map_coords = new Point[num_map_points];
    SignalCovarParams * const s_params =
        new SignalCovarParams(corr_var_s, corr_l_perp, corr_l_para);
    PCGParams * const pcg_params = new PCGParams;

    // Init pixel data.
    puts("Reading pixel data.");
    read_pixel_data(pixel_data_path, num_pixels, pixels);

    puts("Preparing for x solve.");

    // Init map coords.
    init_uniform_map_coords(lx, ly, lz, map_nx, map_ny, map_nz, map_coords);

    // Init PCG Params.
    pcg_params->max_iter = pcg_max_iter;
    pcg_params->step_r = pcg_step_r;
    pcg_params->tol = pcg_tol;

    //
    // Reconstruction
    //

    // convert pixel format.
    WPixel *wpixels = (WPixel *) pixels;
    pixel_to_wpixel(num_pixels, pixels, wpixels);

    // Create x solver.
    WFX_W *wfx = new WFX_W(num_pixels, wpixels, s_params);

    // x solve
    double *x = new double[num_pixels];
    for (int i = 0; i < num_pixels; ++i) { x[i] = 0.0; }

    puts("Starting solve.");

    if (solve_via_pcg) {
        // true for verbose option.
        wfx->solve_pcg(x, pcg_params, true);
    }
    else {
        // true for verbose and exact.
        wfx->solve_cf(x, true, true);
    }

    // map multiply
    // make pixel coords
    const Point * const pixel_coords =
        points_from_pixels_alloc(num_pixels, pixels);

    double *m = new double[num_map_points];
    smp_product(num_pixels, pixel_coords, num_map_points, map_coords,
        s_params, x, m);

    // Write map field.
    printf("Writing map file %s.\n", map_path.c_str());
    FILE *map_file = fopen(map_path.c_str(), "w");
    fwrite(m, sizeof(double), num_map_points, map_file);
    fclose(map_file);

    /*
    if (option_noise_covar) {
        puts("Computing covariance in N >> S limit.");
        //ds_smp_w_spm(map);
        map_file = fopen("map_covar_nd.bin", "w");
        fwrite(m, sizeof(double), num_map_points, map_file);
        fclose(map_file);
    }

    if (option_map_covar) {
        puts("Computing covariance diag.");
        //ds_map_covar_diag(map);
        map_file = fopen("map_covar.bin", "w");
        fwrite(m, sizeof(double), num_map_points, map_file);
        fclose(map_file);
    }
    */

    delete pcg_params;
    delete s_params;
    delete [] map_coords;
    delete [] pixels;

    // report timing.
    printf("Total time: %.2f s.\n", total_timer.elapsed());

    return 0;
}
