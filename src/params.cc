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

    p.num_skewers = 10;
    p.num_pixels = 10;
    p.pix_n = 100;
    p.pix_dz = 1.0;

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

    p.option_noise_covar = 1;
    p.option_map_covar = 0;
    p.option_smooth_map = 0;

    p.skewer_x_path = "skewer_x.bin";
    p.skewer_y_path = "skewer_y.bin";
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
        else if (strcmp(key, "num_skewers") == 0) {
            p.num_skewers = atoi(value);
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
            p.corr_l_para_2 = x*x;
        }
        else if (strcmp(key, "corr_l_perp") == 0) {
            double x = atof(value);
            p.corr_l_perp_2 = x*x;
        }
        else if (strcmp(key, "pcg_tol") == 0) {
            p.pcg_tol = atof(value);
        }
        else if (strcmp(key, "pcg_max_iter") == 0) {
            p.pcg_max_iter = atoi(value);
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
        else if (strcmp(key, "skewer_x_path") == 0) {
            p.skewer_x_path = value;
        }
        else if (strcmp(key, "skewer_y_path") == 0) {
            p.skewer_y_path = value;
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

    p.pix_n = p.num_pixels * p.num_skewers;
    p.pix_dz = p.lz / p.num_pixels;
    p.map_n = p.map_nx * p.map_ny * p.map_nz;
    p.map_dx = p.lx / p.map_nx;
    p.map_dy = p.ly / p.map_ny;
    p.map_dz = p.lz / p.map_nz;
}
