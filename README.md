FLIPS MATLAB Package - By M. R. Sheriff, F. F. Redel, and P. Mohajerin Esfahani, see [1] for a detailed description of FLIPS. In particular, the MATLAB code in this sub-branch is for the problem of Compressed Sensing, see [1, Section 4.2] for more details. We solve the LIP

min_f           ||f||_1

s.t.            ||x – (CD)f||_2 <= epsilon ,

where D is the inverse-DCT matrix and C is the randomly generated measurement matrix.

1. The "Main_Compressed_Sensing.m" is the main file that generates the problem data for a given image and calls the FLIPS solver to solve the corresponding LIP

2. The "The "FLIPS_Solver.m" file contains the FLIPS solver.

3. The "show_results.m" file reconstructs the image from the solution from the FLIPS solver


These files are allowed to be adjusted. However, it is not allowed to publish or distribute these files without permission from the authors.

This research was supported by the European Research Council (ERC) under the European Unions Horizon 2020 research and innovation program (TRUST-949796).

%%%=======================================================================================

-References


[1] M. R. Sheriff, F. F. Redel, and P. Mohajerin Esfahani, Fast Algorithm for Constrained Linear Inverse Problems. arXiv: 2212.01068.


