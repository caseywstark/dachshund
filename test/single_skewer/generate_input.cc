/*
2014, the dachshund authors.
*/

#include <cstdio>

int
main()
{
    int num_skewers = 1;
    int num_pixels = 4;
    int n = num_skewers * num_pixels;

    // set coords
    double x[num_skewers];
    double y[num_skewers];
    x[0] = 0.5;
    y[0] = 0.5;

    double data[4] = {0.0, 1.0, 1.0, 0.0};

    double w[4];
    for (int i = 0; i < n; ++i) {
        w[i] = 1.0;
    }

    FILE *coord_file = fopen("skewer_coords.bin", "w");
    fwrite(x, sizeof(double), num_skewers, coord_file);
    fwrite(y, sizeof(double), num_skewers, coord_file);
    fclose(coord_file);

    FILE *data_file = fopen("skewer_data.bin", "w");
    fwrite(data, sizeof(double), n, data_file);
    fclose(data_file);

    FILE *w_file = fopen("skewer_weights.bin", "w");
    fwrite(w, sizeof(double), n, w_file);
    fclose(w_file);

    return 0;
}
