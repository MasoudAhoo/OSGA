OSGA package
=============


This matlab package provides a generic solver for Optimal SubGradient Algorithm (OSGA) for solving structured nonsmooth convex problems of the form

  min_x sum_{i=1}^{n_1} f_i(A_i x)+sum_{j=1}^{n_2} phi_j(W_j x), x in C

 for proper, convex, and lower semicontinuous f_i (i=1,...,n_1) and g_j (j=1,...,n_2), and C is a simple constraint.

We provide a solver for solving the above-mentioned problem:
- OSGA: Optimal SubGradient Algorithm.

# Using the package

## Installation

All the functions in the package are given in the folder /MatlabCodes.

In order to use the function we recommend to execute the following command

```Matlab
addpath(genpath('.'))
```

if you are not working in the root folder of the package or replacing '.' by the location of the folder on your machine.


## Examples:

We recommend to look at the following files to see how to use the package:
* Demos/Demo_deblurring.m: Image deblurring;
* Demos/Demo_deblurring_nonngativeCon.m: Image deblurring with nonnegativity constraints;
* Demos/Demo_denoising.m: Image denoising;
* Demos/Demo_inpainting.m: Image inpainting;
* Demos/Demo_L1.m: l1-minimization;
* Demos/Demo_L22L22R_BoundCon.m: Bound constrained L22L22R minimization;
* Demos/Demo_Ridge_regression.m: Ridge regression;
* Demos/Demo_sparse_recovery.m: sparse recovery.

## Solving your own optimization problem

please see the user manual: User's_manual_for_OSGA.pdf.

# References

[1] M. Ahookhosh, Optimal subgradient methods: computational properties for
large-scale linear inverse problems. Optimization and Engineering, 19(4), 815-
844 (2018).

[2] M. Ahookhosh, A. Neumaier, An optimal subgradient algorithm for large-scale
bound-constrained convex optimization. Mathematical Methods of Operations
Research, 86(1), 123-147 (2017).

[3] M. Ahookhosh, A. Neumaier, Optimal subgradient algorithms for large-scale
convex optimization in simple domains, Numerical Algorithms, 76(4), 1071-
1097 (2017).

[4] M. Ahookhosh, A. Neumaier, An optimal subgradient algorithm with subspace
search for costly convex optimization problems, Bulletin of the Iranian Mathe-
matical Society, 45(3), 883–910 (2019).

[5] M. Ahookhosh, A. Neumaier, Solving structured nonsmooth convex optimiza-
tion with complexity O(ε −1/2 ), Top, 26(1), 110-145, (2018).

[6] A. Neumaier, OSGA: a fast subgradient algorithm with optimal complexity,
Mathematical Programming, 158(1-2), 1-21 (2016).



