
#include <cmath>
#include <cstdlib>

#include "test_utils.h"

bool equal_within_tol(const double a, const double b, const double abs_tol,
        const double rel_tol) {
    // a is the 'expected'.
    return fabs(a - b) <= (abs_tol + rel_tol * fabs(a));
}

void init_rng() {
    srand(time(0));
}

double rng() {
    return rand() / static_cast<double>(RAND_MAX);
}
