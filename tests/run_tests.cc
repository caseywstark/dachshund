
#include <cstdlib>

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

int main( int argc, char* const argv[] )
{
    // global setup...

    // setup rng.
    srand(time(0));

    int result = Catch::Session().run( argc, argv );

    // global clean-up...

    return result;
}
