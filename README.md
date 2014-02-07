Dachshund
=========

Given skewers (delta_F, noise), reconstruct the 3d delta_F field.

Install
-------

Dachshund is a C++ code. You will need a C++ compiler with openmp support.

The make system is straightforward. The root `Makefile` includes the
`platform.make` file, which controls some system-specific settings. There are
some examples in the platforms directory.

Usage
-----

The Dachshund executable is very simple for now. It takes in a file describing
the skewer coordinates and pixel fluxes and noise. It outputs the binary blob
of the 3d reconstruction.

$ ./dachshund input.cfg
