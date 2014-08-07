Dachshund
=========

Wiener Filter reconstruction of 3D Lyman-alpha Forest flux fields.

Building
--------

Dachshund is written in C++. You will need a C++ compiler. OpenMP support is
preferable for performance, but it's optional.

The make system is straightforward. The root `Makefile` includes the
`platform.make` file, which controls some system-specific settings. There are
some examples in the `platforms` directory.

The default `make` target builds the dachshund library, the dachshund
application, and runs the test suite.

Before building for the first time, please make the UnitTest++ library in
`tests/UnitTest++`. This is not required for the library and application, but it
is required to run the tests.

Usage
-----

The dachshund application requires an input config file specifying the run-time
parameters and an input data file containing the pixel data. The config file
format might change during development, so please check against
`app/dachshund.cc` (the main application source) to make sure it is up to date.
Then run the dachshund application with:

    $ ./dachshund.ex input.cfg

The pixel data file should contain the spatial coordinate $(x, y, z)$, delta_F
value $(d)$, and delta_F noise estimate $(n)$ for each pixel in double
precision. If the data vector for each pixel is $p_i = (x_i, y_i, z_i, n_i,
d_i)$, the data file should be a binary blob of $p_0, p_1, \ldots, p_{n_{\rm
pixel} - 1}$.
