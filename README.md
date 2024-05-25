FLIPS MATLAB Package - By M. R. Sheriff, F. F. Redel, and P. Mohajerin Esfahani, see [1] for a detailed description of FLIPS. In particular, the MATLAB code in this sub-branch is for the image-denoising problem that solves

min_f           ||f||_1

s.t.            ||x – Df||_2 <= epsilon,

where D is the inverse-DCT matrix.

1. The "Main_Image_Denoising.m" file is the main file that extracts patches from the image and calls various algorithms to solve the LIPs.

2. The "FLIPS_denoising_final.m" file is used to find an optimal solution of the LIPs corresponding to every patch.

3. The "show_results.m" file computes the recovered images for all algorithms

6. The "show_results.m" file computes the recovered images from the solution to all the patches and also creates the convergence plots.


These files are allowed to be adjusted. However, it is not allowed to publish or distribute these files without permission from the authors. 

This research was supported by the European Research Council (ERC) under the European Unions Horizon 2020 research and innovation program (TRUST-949796).


%%%==============================================================================

-References


[1] M. R. Sheriff, F. F. Redel, and P. Mohajerin Esfahani, Fast Algorithm for Constrained Linear Inverse Problems. arXiv: 2212.01068.
