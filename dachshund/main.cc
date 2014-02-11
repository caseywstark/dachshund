/*
2014, the dachshund authors.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>

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
    // Fancy command line parsing.
    if (argc != 2) {
        fputs("Bad number of arguments. Expected 1.\n", stderr);
        print_usage();
        exit(1);
    }

    std::string input_path = argv[1];

    //
    // Set default params.
    //

    int num_skewers = 100;
    int num_pixels = 100;
    std::string coords_path = "skewer_coords.bin";
    std::string data_path = "skewer_data.bin";
    std::string weights_path = "skewer_weights.bin";
    int nx = 64;
    int ny = 64;
    int nz = 64;
    double lx = 100.0;
    double ly = 100.0;
    double lz = 100.0;
    std::string output_path = "map.bin";
    double sigma_F = 0.1;
    double l_perp = 0.5;
    double l_para = 1.0;
    double pcg_tol = 1.0e-3;
    int pcg_max_iter = 100;

    //
    // Parse input/config.
    //

    const static int max_line_length = 1000;
    char config_buf[max_line_length];
    char key[100];
    char value[100];

    // Open input file.
    FILE *config_file = fopen(input_path.c_str(), "r");
    // Make sure the input file is ok.
    if (config_file == NULL) {
        fprintf(stderr, "Could not load file %s.\n", input_path.c_str());
        exit(1);
    }

    printf("Reading config file %s.\n", input_path.c_str());

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
            num_skewers = atoi(value);
        }
        else if (strcmp(key, "num_pixels") == 0) {
            num_pixels = atoi(value);
        }
        else if (strcmp(key, "coords_path") == 0) {
            coords_path = value;
        }
        else if (strcmp(key, "data_path") == 0) {
            data_path = value;
        }
        else if (strcmp(key, "weights_path") == 0) {
            weights_path = value;
        }
        else if (strcmp(key, "nx") == 0) {
            nx = atoi(value);
        }
        else if (strcmp(key, "ny") == 0) {
            ny = atoi(value);
        }
        else if (strcmp(key, "nz") == 0) {
            nz = atoi(value);
        }
        else if (strcmp(key, "lx") == 0) {
            lx = atof(value);
        }
        else if (strcmp(key, "ly") == 0) {
            ly = atof(value);
        }
        else if (strcmp(key, "lz") == 0) {
            lz = atof(value);
        }
        else if (strcmp(key, "output_path") == 0) {
            output_path = value;
        }
        else if (strcmp(key, "sigma_F") == 0) {
            sigma_F = atof(value);
        }
        else if (strcmp(key, "l_perp") == 0) {
            l_perp = atof(value);
        }
        else if (strcmp(key, "l_para") == 0) {
            l_para = atof(value);
        }
        else if (strcmp(key, "pcg_tol") == 0) {
            pcg_tol = atof(value);
        }
        else if (strcmp(key, "pcg_max_iter") == 0) {
            pcg_max_iter = atoi(value);
        }
        else {
            fprintf(stderr, "Unknown key '%s' in config.\n", key);
            exit(1);
        }
    }

    int64_t num_cells = (int64_t)nx * (int64_t)ny * (int64_t)nz;
    int64_t num_pixel_elements = (int64_t)num_pixels * (int64_t)num_skewers;
    double dx = lx / nx;
    double dy = ly / ny;
    double dz = lz / nz;
    double dz_pix = lz / num_pixels;

    double *skewer_x = new double[num_skewers];
    double *skewer_y = new double[num_skewers];
    double *data_vec = new double[num_pixel_elements];
    double *weights_vec = new double[num_pixel_elements];

    //
    // Read in skewers.
    //

    puts("Reading skewer files.");

    FILE *coords_file = fopen(coords_path.c_str(), "r");
    if (coords_file == NULL) {
        fprintf(stderr, "Could not load file %s.\n", coords_path.c_str());
        exit(1);
    }

    fread(skewer_x, sizeof(double), num_skewers, coords_file);
    fread(skewer_y, sizeof(double), num_skewers, coords_file);
    fclose(coords_file);

    FILE *data_file = fopen(data_path.c_str(), "r");
    if (data_file == NULL) {
        fprintf(stderr, "Could not load file %s.\n", data_path.c_str());
        exit(1);
    }
    fread(data_vec, sizeof(double), num_pixel_elements, data_file);
    fclose(data_file);

    FILE *weights_file = fopen(weights_path.c_str(), "r");
    if (weights_file == NULL) {
        fprintf(stderr, "Could not load file %s.\n", weights_path.c_str());
        exit(1);
    }
    fread(weights_vec, sizeof(double), num_pixel_elements, weights_file);
    fclose(weights_file);

    // temporary
    double *noise_vec = weights_vec;
    for (int i = 0; i < num_pixel_elements; ++i) {
        noise_vec[i] = 1.0 / noise_vec[i];
    }

    // DEBUG
    // Check input data
    check_finite(num_skewers, skewer_x);
    check_finite(num_skewers, skewer_y);
    check_finite(num_skewers, data_vec);
    check_finite(num_skewers, weights_vec);

    // Create gpt.
    int gt_n = 100;
    double gt_dx = 4.0 / gt_n;
    GT *gt = new GT(gt_n, gt_dx);

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
    p1->gt = gt;

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
    p2->gt = gt;

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
