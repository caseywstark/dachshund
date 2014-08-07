Dachshund design notes
======================

These are some notes for other developers. Apologies if these are not kept up
to date...

Library organization
--------------------

There are two main steps to computing the reconstructed map $m = S x = S (S +
N)^{-1} d$. The first is the solve for $x$, and the second is the $S x$ product.
Most of the code is written for the first step.

...TODO...

Misc.
-----

- When using MPI, each rank has a full copy of the data. The parallelization
is only splitting the matrix-vector multiplies, as that is the dominant
operation in each step.

