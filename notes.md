Source notes
============

Some notes for other developers.

Layout:
- `app` contains application sources (the compilation unit with a `main`). Currently just the main dachshund app.
- `eigen` third-party source.
- `examples` contains an example problem.
- `lib` is the main source that makes up the library.
- `tests` contains third-party Catch source, the test suite, and the test runner.

There are two main steps to computing the reconstructed map $m = S x = S (S +
N)^{-1} d$. The first is the solve for $x$, and the second is the $S x$ product.
Most of the code is devoted to the first step, contained in `map.cc`.

If we need to run larger problems in the future, we can extend to MPI
parallelism. Each MPI rank can store of a full copy of the data.
The parallelization is only splitting the matrix-vector multiplies, where each
rank would be responsible for a chunk of rows and then combine their results.
