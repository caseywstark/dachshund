Source design
=============

Some notes for other developers.

There are two main steps to computing the reconstructed map $m = S x = S (S +
N)^{-1} d$. The first is the solve for $x$, and the second is the $S x$ product.
Most of the code is devoted to the first step.



If we need to run larger problems in the future, we can extend to MPI
parallelism. Each MPI rank can store of a full copy of the data.
The parallelization is only splitting the matrix-vector multiplies, where each
rank would be responsible for a chunk of rows and then combine their results.

