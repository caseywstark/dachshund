/*
2014, the dachshund authors.

Very simple module -- just declare the parameter object so we can access it
throughout the code.
*/

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "params.h"

DSParams p;

void
ds_params_init(const std::string config_path)
{
    // Default values.
    p.lx = 10.0;
    p.ly = 10.0;
    p.lz = 10.0;

    p.num_pixels = 100;

    p.map_nx = 10;
    p.map_ny = 10;
    p.map_nz = 10;
    p.map_n = 1000;
    p.map_dx = 1.0;
    p.map_dy = 1.0;
    p.map_dz = 1.0;

    p.corr_var_s = 1.0;
    p.corr_l_para_2 = 1.0;
    p.corr_l_perp_2 = 1.0;

    p.pcg_max_iter = 1;
    p.pcg_tol = 1.0;
    p.pcg_step_r = 10;

    p.option_noise_covar = 0;
    p.option_map_covar = 0;
    p.option_smooth_map = 0;

    p.pixel_x_path = "pixel_x.bin";
    p.pixel_y_path = "pixel_y.bin";
    p.pixel_z_path = "pixel_z.bin";
    p.pixel_data_path = "pixel_data.bin";
    p.pixel_weights_path = "pixel_weights.bin";
    p.map_path = "map.bin";

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
            p.lx = atof(value);
        }
        else if (strcmp(key, "ly") == 0) {
            p.ly = atof(value);
        }
        else if (strcmp(key, "lz") == 0) {
            p.lz = atof(value);
        }
        else if (strcmp(key, "num_pixels") == 0) {
            p.num_pixels = atoi(value);
        }
        else if (strcmp(key, "map_nx") == 0) {
            p.map_nx = atoi(value);
        }
        else if (strcmp(key, "map_ny") == 0) {
            p.map_ny = atoi(value);
        }
        else if (strcmp(key, "map_nz") == 0) {
            p.map_nz = atoi(value);
        }
        else if (strcmp(key, "corr_var_s") == 0) {
            p.corr_var_s = atof(value);
        }
        else if (strcmp(key, "corr_l_para") == 0) {
            double x = atof(value);
            p.corr_l_para_2 = 2.0*x*x;
        }
        else if (strcmp(key, "corr_l_perp") == 0) {
            double x = atof(value);
            p.corr_l_perp_2 = 2.0*x*x;
        }
        else if (strcmp(key, "pcg_tol") == 0) {
            p.pcg_tol = atof(value);
        }
        else if (strcmp(key, "pcg_max_iter") == 0) {
            p.pcg_max_iter = atoi(value);
        }
        else if (strcmp(key, "pcg_step_r") == 0) {
            p.pcg_step_r = atoi(value);
        }
        else if (strcmp(key, "option_noise_covar") == 0) {
            p.option_noise_covar = atoi(value);
        }
        else if (strcmp(key, "option_map_covar") == 0) {
            p.option_map_covar = atoi(value);
        }
        else if (strcmp(key, "option_smooth_map") == 0) {
            p.option_smooth_map = atoi(value);
        }
        else if (strcmp(key, "pixel_x_path") == 0) {
            p.pixel_x_path = value;
        }
        else if (strcmp(key, "pixel_y_path") == 0) {
            p.pixel_y_path = value;
        }
        else if (strcmp(key, "pixel_z_path") == 0) {
            p.pixel_z_path = value;
        }
        else if (strcmp(key, "pixel_data_path") == 0) {
            p.pixel_data_path = value;
        }
        else if (strcmp(key, "pixel_weights_path") == 0) {
            p.pixel_weights_path = value;
        }
        else if (strcmp(key, "map_path") == 0) {
            p.map_path = value;
        }
        else {
            fprintf(stderr, "Unknown key '%s' in config.\n", key);
            exit(1);
        }
    }

    fclose(config_file);

    p.map_n = p.map_nx * p.map_ny * p.map_nz;
    p.map_dx = p.lx / p.map_nx;
    p.map_dy = p.ly / p.map_ny;
    p.map_dz = p.lz / p.map_nz;
}

void
ds_params_print()
{
    printf("\n");
    printf("Running dachshund with:\n");
    printf("  pixels/map:\n");
    printf("    lx ly lz = %f %f %f\n", p.lx, p.ly, p.lz);
    printf("    npix = %i\n", p.num_pixels);
    printf("    map nx ny nz = %i %i %i, n = %i\n", p.map_nx, p.map_ny, p.map_nz, p.map_n);
    printf("    map dx dy dz = %f %f %f\n", p.map_dx, p.map_dy, p.map_dz);
    printf("  signal:\n");
    printf("    corr var_s = %f\n", p.corr_var_s);
    printf("    2*(l_perp)^2 = %f    2*(l_para)^2 = %f\n", p.corr_l_perp_2, p.corr_l_para_2);
    printf("  PCG:\n");
    printf("    max iter = %i    step r = %i\n", p.pcg_max_iter, p.pcg_step_r);
    printf("    tol = %f\n", p.pcg_tol);
    printf("  options:\n");
    printf("    map covar = %i\n", p.option_map_covar);
    printf("    noise covar = %i\n", p.option_noise_covar);
    printf("    smooth map = %i\n", p.option_smooth_map);
    printf("  paths:\n");
    printf("    x = %s\n", p.pixel_x_path.c_str());
    printf("    y = %s\n", p.pixel_y_path.c_str());
    printf("    z = %s\n", p.pixel_z_path.c_str());
    printf("    data = %s\n", p.pixel_data_path.c_str());
    printf("    weights = %s\n", p.pixel_weights_path.c_str());
    printf("    map = %s\n", p.map_path.c_str());
    printf("\n");
}
