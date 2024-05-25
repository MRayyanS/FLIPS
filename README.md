FLIPS MATLAB Package - By M. R. Sheriff, F. F. Redel, and P. Mohajerin Esfahani, see [1] for a detailed description of FLIPS. In particular, the MATLAB code in this sub-branch is for the image-denoising problem that solves

min_f           ||f||_1

s.t.            ||x – Df||_2 <= epsilon,

where D is the inverse-DCT matrix. We also compare FLIPS with other standard algorithms namely C-SALSA and Chambolle-Pock algorithms.

1. The "Main_Image_Denoising_CPUtimes.m" file is the main file that extracts patches from the image and calls various algorithms to solve the LIPs.

2. The "FLIPS_denoising_final.m" file is used to find an optimal solution of the LIPs corresponding to each patch which is used for comparison with other algorithms.

3. The "FLIPS_CPU_time.m" file solves the LIPs using the FLIPS algorithm [1] and records the corresponding CPU times

4. The "C_SALSA_CPU_time.m" file solves the LIPs using the C-SALSA algorithm [2] and records the corresponding CPU times

5. The "ChambollePock_CPU_time.m" file solves the LIPs using the Chambolle-Pock algorithm [3] and records the corresponding CPU times

6. The "show_results.m" file computes the recovered images for all algorithms


These files are allowed to be adjusted. However, it is not allowed to publish or distribute these files without permission from the authors. 

This research was supported by the European Research Council (ERC) under the European Unions Horizon 2020 research and innovation program (TRUST-949796).


%%%=================================================================
-References
%%%=================================================================


[1] M. R. Sheriff, F. F. Redel, and P. Mohajerin Esfahani, Fast Algorithm for Constrained Linear Inverse Problems. arXiv: 2212.01068.

[2] M. V. Afonso, J. M. Bioucas-Dias, and M. A. T. Figueiredo, An Augmented Lagrangian Approach to the Constrained Optimization Formulation of Imaging Inverse Problems, IEEE Transactions On Image Processing, (2009).
    
[3] A. Chambolle and T. Pock. On the Ergodic Convergence rates of a First-order Primal-Dual Algorithm. Mathematical Programming, 159(1-2):253–287, 9 2016




