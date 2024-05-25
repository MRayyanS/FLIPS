FLIPS MATLAB Package - By M. R. Sheriff, F. F. Redel, and P. Mohajerin Esfahani. This is a MATLAB repository for the Fast Linear Inverse Problems Solver (FLIPS) introduced in [1].

FLIPS is an algorithm to solve the constrained linear inverse problems (LIP), see [1] for a detailed description of FLIPS. In particular, the software here is a MATLAB package for the image-denoising problem that solves

min_f           c(f)

s.t.            ||x – Df||_2 <= epsilon

The main branch "main_FLIPS_solver" contains a generic FLIPS solver module and sample images to test examples of linear inverse problems. Each of the sub-branches corresponds to a specific instance of LI with different choices of the cost function c() and linear map D. All the files needed to run for that example of LIP are contained in that sub-branch except the example figures. Please read the readme files of the respective sub-branch for more explanation.

These files are allowed to be adjusted. However, without the author's permission, it is not allowed to publish or distribute these files. 

This research was supported by the European Research Council (ERC) under the European Unions Horizon 2020 research and innovation program (TRUST-949796).



%%%=======================================================================================

-  References


[1] M. R. Sheriff, F. F. Redel, and P. Mohajerin Esfahani, Fast Algorithm for Constrained Linear Inverse Problems. arXiv: 2212.01068.



