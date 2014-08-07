#ifndef __pm_test_util_h__
#define __pm_test_util_h__

#include <ctime>
#include <cstdlib>

#include <UnitTest++.h>

#define CHECK_DTOL(a, b, atol, rtol) CHECK(check_double_tol(a, b, atol, rtol))
#define CHECK_ARRAY_DTOL(n, a, b, atol, rtol) \
    CHECK(check_array_double_tol(n, a, b, atol, rtol))

bool
check_double_tol(const double a, const double b, const double atol,
    const double rtol);

bool
check_array_double_tol(const int n, const double * const a,
    const double * const b, const double atol, const double rtol);

double
rng();

#endif
