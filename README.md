Dachshund
=========

Tomographic reconstruction of Lyman-alpha forest flux using a Wiener filter.

Building
--------

Dachshund is written in C++. You will need a C++ compiler, preferably with
OpenMP support.

The make system is straightforward. The root `Makefile` includes the
`platform.make` file, which controls some system-specific settings. There are
some examples in the `platforms` directory.

The default `make` target builds the dachshund library, the dachshund
application, and runs the test suite.

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

Dependencies
------------

Dachshund uses [Eigen](http://eigen.tuxfamily.org/) for some linear algebra
routines and [Catch](https://github.com/philsquared/Catch) for organizing tests.
Both are included in the repository, so there is no need to get them separately.
