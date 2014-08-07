/*
The test runner main.
*/

#include <ctime>
#include <cstdlib>

#include <UnitTest++.h>

int main()
{
    // setup rng.
    srand(time(0));

    // will eventually do mpi stuff here...
    return UnitTest::RunAllTests();
}
