
#include <math.h>

#include "test_utils.h"

bool
check_double_tol(const double a, const double b, const double atol,
    const double rtol)
{
    return fabs(a - b) <= (atol + rtol * fabs(b));
}

bool
check_array_double_tol(const int n, const double * const a,
    const double * const b, const double atol, const double rtol)
{
    bool close = true;
    for (int i = 0; i < n; ++i) {
        close &= check_double_tol(a[i], b[i], atol, rtol);
    }
    return close;
}

double
rng()
{
    return rand() / double(RAND_MAX);
}
