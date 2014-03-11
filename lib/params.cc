/*
2014, the dachshund authors.

Very simple module -- just declare the parameter object so we can access it
throughout the code.
*/

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "params.h"

// Init params object.
// Unfortunately, we need static initialization here... getting messy...
DSParams p = {
    10, 10, 10, 1000, // nx, ny, nz, num_cells
    10.0, 10.0, 10.0, //lx, ly, lz
    1.0, 1.0, 1.0, //dx, dy, dz
    100, //num_skewers
    10, //num_pixels
    1000, //num_pixel_elements
    1.0, //dz_pix
    1.0, //l_perp
    1.0, // l_para
    1.0, // sigma
    10,  //pcg_max_iter
    1.0, // pcg_tol
    0, //compute_covar
    NULL, //pixels
    NULL, //map
    "skewer_x.bin", // skewer_x_path;
    "skewer_y.bin", //skewer_y_path;
    "pixel_data.bin", //pixel_data_path;
    "pixel_weights.bin", //pixel_weights_path;
    "map.bin", //map_path;
    NULL //gt;
};

void
ds_params_init(const std::string config_path)
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

        if (strcmp(key, "num_skewers") == 0) {
            p.num_skewers = atoi(value);
        }
        else if (strcmp(key, "num_pixels") == 0) {
            p.num_pixels = atoi(value);
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
        else if (strcmp(key, "nx") == 0) {
            p.nx = atoi(value);
        }
        else if (strcmp(key, "ny") == 0) {
            p.ny = atoi(value);
        }
        else if (strcmp(key, "nz") == 0) {
            p.nz = atoi(value);
        }
        else if (strcmp(key, "lx") == 0) {
            p.lx = atof(value);
        }
        else if (strcmp(key, "ly") == 0) {
            p.ly = atof(value);
        }
        else if (strcmp(key, "lz") == 0) {
            p.lz = atof(value);
        }
        else if (strcmp(key, "sigma_F") == 0) {
            p.sigma = atof(value);
        }
        else if (strcmp(key, "l_perp") == 0) {
            p.l_perp = atof(value);
        }
        else if (strcmp(key, "l_para") == 0) {
            p.l_para = atof(value);
        }
        else if (strcmp(key, "pcg_tol") == 0) {
            p.pcg_tol = atof(value);
        }
        else if (strcmp(key, "pcg_max_iter") == 0) {
            p.pcg_max_iter = atoi(value);
        }
        else if (strcmp(key, "compute_covar") == 0) {
            p.compute_covar = atoi(value);
        }
        else {
            fprintf(stderr, "Unknown key '%s' in config.\n", key);
            exit(1);
        }
    }

    fclose(config_file);

    p.num_cells = (int64_t)p.nx * (int64_t)p.ny * (int64_t)p.nz;
    p.num_pixel_elements = (int64_t)p.num_pixels * (int64_t)p.num_skewers;
    p.dx = p.lx / p.nx;
    p.dy = p.ly / p.ny;
    p.dz = p.lz / p.nz;
    p.dz_pix = p.lz / p.num_pixels;
}
