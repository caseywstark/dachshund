/*
2014, the dachshund authors.
*/

#ifndef __DS_LUM_H__
#define __DS_LUM_H__

#include "gt.h"

class LUM {
  public:
    LUM(int n, int m, double (*element_func)(int i, int j, void *params),
        void *params);

    // members
    int n;
    int m;
    double (*element_func)(int i, int j, void *params);
    void *params;

    // methods
    double operator()(int i, int j);

  private:
    //
};

struct smd_params_struct {
    int nx;
    int ny;
    int nz;
    double dx;
    double dy;
    double dz;

    int ns;
    int np;
    double dz_pix;
    double *skewer_x;
    double *skewer_y;

    double l_perp;
    double l_para;
    double sigma;

    GT *gt;
};

typedef struct smd_params_struct SMDParams;

double
smd_element_func(int i, int j, void *params);

struct sddn_params_struct {
    int ns;
    int np;
    double dz_pix;
    double *skewer_x;
    double *skewer_y;

    double l_perp;
    double l_para;
    double sigma;

    double *n;
    GT *gt;
};

typedef struct sddn_params_struct SDDNParams;

double
sddn_element_func(int i, int j, void *params);

#endif
