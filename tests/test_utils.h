#ifndef __ds_test_utils_h__
#define __ds_test_utils_h__

#include "gtest/gtest.h"

#define ASSERT_EQTOL(expected, actual, abs_tol, rel_tol) \
    ASSERT_PRED4(equal_within_tol, expected, actual, abs_tol, rel_tol)
#define EXPECT_EQTOL(expected, actual, abs_tol, rel_tol) \
    EXPECT_PRED4(equal_within_tol, expected, actual, abs_tol, rel_tol)

bool equal_within_tol(const double a, const double b, const double abs_tol,
    const double rel_tol);

void init_rng();

double rng();


#endif
