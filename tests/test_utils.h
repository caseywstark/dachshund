#ifndef __ds_test_utils_h__
#define __ds_test_utils_h__

#include <cmath>
#include <cstdlib>

static inline bool
dtol(const double a, const double b, const double atol,
    const double rtol)
{
    return fabs(a - b) <= (atol + rtol * fabs(b));
}

static inline double
rng()
{
    return rand() / double(RAND_MAX);
}

#endif
