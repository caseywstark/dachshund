/*
2014, the dachshund authors.
*/

#include <math.h>
#include <cstdio>

#include "util.h"

void
check_finite(const int n, const double * const v)
{
    for (int i = 0; i < n; ++i) {
        bool f = isfinite(v[i]);
        if (!f) {
            printf("non-finite value v[%i] = %e\n", i, v[i]);
        }
    }
}
