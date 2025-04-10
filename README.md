# nearest-singular-matrix-valued-function
This repository contains the code for the numerical approximation of the distance to singularity for matrix-valued functions, as described in

M. Gnazzo, N. Guglielmi. "On the numerical approximation of the distance to singularity for matrix-valued functions", available on [arXiv](https://arxiv.org/abs/2309.01220).

Given a regular matrix-valued function in the form $\sum_{i=1}^k f_i(\lambda) A_i$, we approximate closest singular matrix-valued function $\sum_{i=1}^k f_i(\lambda) (A_i + \Delta A_i)$, perturbing the coefficients $A_i$, for $i=1,\ldots,k$. It can be applied both to matrix polynomials and matrix-valued functions, with nonlinearities different from the polynomial one.

## How to use it:
* **struct_distance**: contains the function for the numerical approximation of the distance to singularity; it needs a cell array with the coeffiecients $A_i$ and a function handle for $f_i$.
* **Script_delay** and **Script_polynomial**: contain two examples and can be used as a reference to run the method.
