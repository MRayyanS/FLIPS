FLIPS MATLAB Package - By M. R. Sheriff, F. F. Redel, and P. Mohajerin Esfahani, see [1] for a detailed description of FLIPS. In particular, the MATLAB code in this sub-branch is for the Binary Selection problem that solves

min_f           ||f||_inf

s.t.            ||x – Df||_2 <= epsilon .

For a better understanding of the Binary selection problem, see Section 4.1 of [1].

1. The "Main_Binary_Selection.m" file is the main file that generates synthetic problem data for the LIP and then calls FLIPS to solve it.
 
2. The "FLIPS_Solver.m" file contains the FLIPS solver

3. The "show_results.m" creates the convergence plots and converts them into pdfs


These files are allowed to be adjusted. However, it is not allowed to publish or distribute these files without permission from the authors.

This research was supported by the European Research Council (ERC) under the European Unions Horizon 2020 research and innovation program (TRUST-949796).



%%%=======================================================================================

-References


[1] M. R. Sheriff, F. F. Redel, and P. Mohajerin Esfahani, Fast Algorithm for Constrained Linear Inverse Problems. arXiv: 2212.01068.
