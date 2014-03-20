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
    check(skewer_file, "Could not open file %s.", p.skewer_x_path.c_str());
    fread(skewer_x, sizeof(double), p.num_skewers, skewer_file);
    fclose(skewer_file);

    skewer_file = fopen(p.skewer_y_path.c_str(), "r");
    check(skewer_file, "Could not open file %s.", p.skewer_y_path.c_str());
    fread(skewer_y, sizeof(double), p.num_skewers, skewer_file);
    fclose(skewer_file);

    FILE *data_file = fopen(p.pixel_data_path.c_str(), "r");
    check(data_file, "Could not open file %s.", p.pixel_data_path.c_str());
    fread(pixel_data, sizeof(double), p.pix_n, data_file);
    fclose(data_file);

    FILE *weights_file = fopen(p.pixel_weights_path.c_str(), "r");
    check(weights_file, "Could not open file %s.", p.pixel_weights_path.c_str());
    fread(pixel_weights, sizeof(double), p.pix_n, weights_file);
    fclose(weights_file);

error:
    return -1;
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

    pcg_wsppi(p.pix_n, x, b, p.pcg_max_iter, p.pcg_tol, true);

    printf("PCG time: %g ms.\n", pcg_timer->elapsed());

    //
    // Matrix multiply for map.
    //

    // Allocate map.
    double *map = new double[p.map_n];

    printf("Multiply for map_n = %i, total n = %.2e.\n", p.map_n, (double)p.map_n * p.pix_n);

    // m = S^mp x
    // Each thread gets a block of cells.
    ds_smp_x(x, map);

    // Write map field.
    printf("Writing map file %s.\n", p.map_path.c_str());
    FILE *map_file = fopen(p.map_path.c_str(), "w");
    fwrite(map, sizeof(double), p.map_n, map_file);
    fclose(map_file);
/*
    if (p.option_compute_covar) {
        puts("Computing covariance diag.");
        for (int i = 0; i < p.map_n; ++i) {
            printf("cell %i -- %f.\n", i, (double)i/p.map_n);
            // Compute x_i = (S +N)^{-1}_{ij} S^{pm}_{j beta}
            pcg_wsppi_spm(p.pix_n, i, x, p.pcg_max_iter, p.pcg_tol, true);

            // Next step C_{alpha i}  = S^{mp}_{alpha i} x_i
            // Store in the map vector.
            map[i] = 0.0;

            double xi = map_coords[i].x;
            double yi = map_coords[i].y;
            double zi = map_coords[i].z;
            for (j = 0; j < p.pix_n; ++j) {
                // Separation vector components.
                double dx = xi - pixels[j].x;
                double dy = yi - pixels[j].y;
                double dz = zi - pixels[j].z;

                // The gaussian terms.
                double x_perp_2 = (dx*dx + dy*dy) / p.corr_l_perp_2;
                double x_para_2 = dz*dz / p.corr_l_para_2;

                // S = sigma^2 exp(...) exp(...)
                double g_perp = ds_linterp(x_perp_2, gt.n, gt.dx, gt.table);
                double g_para = ds_linterp(x_para_2, gt.n, gt.dx, gt.table);
                double smp_ij = p.corr_var_s * g_perp * g_para;

                map[i] += smp_ij * x[j];
            }
        }

        map_file = fopen("map_covar.bin", "w");
        fwrite(map, sizeof(double), p.map_n, map_file);
        fclose(map_file);
    }
    */

    printf("Total time: %g ms.\n", total_timer->elapsed());

    return 0;
}
