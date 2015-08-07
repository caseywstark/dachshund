/*
The test runner main.
*/

// Cheat to turn gtest into a header-only lib.
#include "gtest/src/gtest-all.cc"
#include "test_utils.h"

int main(int argc, char **argv)
{
    init_rng();
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    return result;
}
