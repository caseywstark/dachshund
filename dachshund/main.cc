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
    Timer *total_timer = new Timer();

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

    // p is the params object declared in the lib.
    p.num_skewers = 100;
    p.num_pixels = 100;
    p.nx = 64;
    p.ny = 64;
    p.nz = 64;
    p.sigma = 1.0;
    p.l_perp = 1.0;
    p.l_para = 1.0;
    p.pcg_tol = 1.0e-3;
    p.pcg_max_iter = 100;

    std::string coords_path = "skewer_coords.bin";
    std::string data_path = "skewer_data.bin";
    std::string weights_path = "skewer_weights.bin";
    std::string output_path = "map.bin";

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
            p.num_skewers = atoi(value);
        }
        else if (strcmp(key, "num_pixels") == 0) {
            p.num_pixels = atoi(value);
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
        else if (strcmp(key, "output_path") == 0) {
            output_path = value;
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
        else {
            fprintf(stderr, "Unknown key '%s' in config.\n", key);
            exit(1);
        }
    }

    fclose(config_file);

    int64_t num_cells = (int64_t)p.nx * (int64_t)p.ny * (int64_t)p.nz;
    int64_t num_pixel_elements = (int64_t)p.num_pixels * (int64_t)p.num_skewers;
    p.dx = p.lx / p.nx;
    p.dy = p.ly / p.ny;
    p.dz = p.lz / p.nz;
    p.dz_pix = p.lz / p.num_pixels;

    double *skewer_x = new double[p.num_skewers];
    double *skewer_y = new double[p.num_skewers];
    double *data_vec = new double[num_pixel_elements];
    double *weights_vec = new double[num_pixel_elements];

    p.skewer_x = skewer_x;
    p.skewer_y = skewer_y;

    //
    // Read in skewers.
    //

    puts("Reading skewer files.");

    FILE *coords_file = fopen(coords_path.c_str(), "r");
    if (coords_file == NULL) {
        fprintf(stderr, "Could not load file %s.\n", coords_path.c_str());
        exit(1);
    }

    fread(skewer_x, sizeof(double), p.num_skewers, coords_file);
    fread(skewer_y, sizeof(double), p.num_skewers, coords_file);
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

    p.noise = noise_vec;

    //
    // Main loop to compute (S + N)^-1 d
    //

    // Create gpt.
    int gt_n = 40;
    double gt_dx = 4.0 / gt_n;
    GT *gt = new GT(gt_n, gt_dx);
    p.gt = gt;

    printf("Entering main loop with n_i = %i, n_j = %i.\n",
        (int)num_cells, (int)num_pixel_elements);

    // We store (S + N)^-1 d here...
    double *b_j = new double[num_pixel_elements];

    // We care about timing starting here.
    Timer *loop_timer = new Timer();

    // Begin parallel block where we give each thread a block of pixels.
    #pragma omp parallel
    {
        int num_threads = omp_get_num_threads();
        printf("Running with %i threads.\n", num_threads);

        int k;
        int n = num_pixel_elements;

        // Allocate thread local arrays.
        double *a_jk = new double[n];
        double *d = new double[n];
        double *r = new double[n];
        double *s = new double[n];
        double *q = new double[n];

        // Set initial a_jk guess.
        for (int k = 0; k < n; ++k) {
            if (k == 0) {
                a_jk[k] = 1.0 / sddn_element_func(k, k);
            }
            else {
                a_jk[k] = 0.0;
            }
        }

        // The main loop over pixels.
        # pragma omp for
        for (int j = 0; j < n; ++j) {

            // DEBUG
            printf("pixel %i\n", (int)j);

            // Setup the residual.
            // r = b - A x

            int jj, kk = 0;
            for (jj = 0; jj < n; ++jj) {
                // Handle b[i]
                if (jj == j) {
                    r[jj] = 1.0;
                }
                else {
                    r[jj] = 0.0;
                }

                //r[i] -= A(i, j) * x[j];
                for (kk = 0; kk < n; ++kk) {
                    r[jj] -= sddn_element_func(jj, kk) * a_jk[kk];
                }
            }

            // The preconditioning.
            // For now, we will just try a Jacobian M (A diag).
            // s = M^-1 r
            for (jj = 0; jj < n; ++jj) {
                d[jj] = r[jj] / sddn_element_func(jj, jj);
            }

            // delta = r^T d
            double delta_new = 0.0;
            for (jj = 0; jj < n; ++jj) {
                delta_new += r[jj] * d[jj];
            }

            double delta_0 = delta_new;

            // Handle singular case.
            //if (delta_0 == 0.0) { return; }

            // the iteration count.
            int iter = 0;
            // the CG loop...
            while (iter < p.pcg_max_iter && delta_new > p.pcg_tol * delta_0) {
                // q = A d
                for (jj = 0; jj < n; ++jj) {
                    q[jj] = 0.0;
                    for (kk = 0; kk < n; ++kk) {
                        q[jj] += sddn_element_func(jj, kk) * d[kk];
                    }
                }

                // alpha = delta_new / (d^T q)
                double denom = 0.0;
                for (jj = 0; jj < n; ++jj) {
                    denom += d[jj] * q[jj];
                }
                double alpha = delta_new / denom;

                // x = x + alpha d
                for (jj = 0; jj < n; ++jj) {
                    a_jk[jj] += alpha * d[jj];
                }

                // Update residual.
                // Approximate update.
                // r = r - alpha q
                for (jj = 0; jj < n; ++jj) {
                    r[jj] -= alpha * q[jj];
                }

                // reapply preconditioner.
                for (jj = 0; jj < n; ++jj) {
                    s[jj] = r[jj] / sddn_element_func(jj, jj);
                }

                // save current delta.
                double delta_old = delta_new;

                // Update delta.
                // delta_new = r^T s
                delta_new = 0.0;
                for (jj = 0; jj < n; ++jj) {
                    delta_new += r[jj] * s[jj];
                }

                // Update d.
                double beta = delta_new / delta_old;
                for (jj = 0; jj < n; ++jj) {
                    d[jj] = s[jj] + beta * d[jj];
                }

                // Finally, update the count.
                ++iter;
            }

            // Dot with data vector.
            double b = 0.0;
            for (k = 0; k < n; ++k) {
                b += a_jk[k] * data_vec[k];
            }

            // Save this pixel result.
            b_j[j] = b;
        } // end omp for

        // don't forget to free the thread work arrays.
        delete [] d;
        delete [] r;
        delete [] q;
        delete [] s;
        delete [] a_jk;

    } // end omp

    printf("Main loop time: %g ms.\n", loop_timer->elapsed());

    delete [] noise_vec;

    //
    // Second loop to compute S b
    //

    // Allocate map.
    double *map_vec = new double[num_cells];

    // Zero out map.
    for (int i = 0; i < num_cells; ++i) {
        map_vec[i] = 0.0;
    }

    // Each thread gets a block of cells.
    #pragma omp parallel for
    for (int i = 0; i < num_cells; ++i) {
        for (int j = 0; j < num_pixel_elements; ++j) {
            map_vec[i] += smd_element_func(i, j) * b_j[j];
        }
    }

    delete [] skewer_x;
    delete [] skewer_y;
    delete [] data_vec;
    delete gt;

    // Write map field.
    printf("Writing map file %s.\n", output_path.c_str());
    FILE *outfile = fopen(output_path.c_str(), "wb");
    fwrite(map_vec, sizeof(double), num_cells, outfile);
    fclose(outfile);

    delete [] map_vec;

    printf("Total time: %g ms.\n", total_timer->elapsed());

    return 0;
}
