FLIPS MATLAB Package - By M. R. Sheriff, F. F. Redel, and P. Mohajerin Esfahani, see [1] for a detailed description of FLIPS. In particular, the MATLAB code in this sub-branch is for the Binary Selection problem that solves

min_f           ||f||_inf

s.t.            ||x – Df||_2 <= epsilon .

For a better understanding of the Binary selection problem, see Section 4.1 of [1].

1. The "Main_Binary_Selection.m" file is the main file that generates synthetic problem data for the LIP and then calls FLIPS to solve the LIPs.

2. The "The "Main_Image_Denoising_CPUtimes.m" file is the main file that extracts patches from the image and calls various algorithms to solve the LIPs."


'Main' file contains the general problem of this package.

The other files contain the following:
- C_SALSA.m --> Contains the C-SALSA solver as described in, see [2] for details of the algorithm 
	
- ChambollePock.m --> Contains the Chambolle-Pock solver, see [3] for details of the algorithm
  
- DCT.m --> The file in which the square DCT-dictionary is computed. 
	
- etafunc.m --> Contains the computations of the cost function eta, as described in [1]
  
- Fista Package --> Package from Tiep Vu https://github.com/tiepvupsu/FISTA open-source (thv102@psu.edu, 4/6/2016).
	
- FLIPS_Solver --> Solver FLIPS for the Quadratic Oracle
	
- Frank_Wolf --> Solver FLIPS for the Linear Oracle
	
- g_descendireciton_FW --> Part of the FLIPS Solver for Linear Oracle.
	
- gradient_eta --> Contains the computations for the gradient of eta, as described [1]
		
- h_updatestep.m --> Update step of the variable h.
	
- inputs --> Contains some standart images that can be used as input. The references of the inputs are listed in the report.
	
- Main.m --> As described above.
	
- patch2image.m --> Function that recreates the image from sliding image patches.
	
- PGD_Oracle.m --> Solver for only Projected Gradient Descent.
 	
- ProjectOntoL1Ball.m --> Projection function to ||.||_1 norm from [4]
	
- soft.m --> Soft thresholding function
	
- stepsize_selection.m --> Exact line search function.

These files are allowed to be adjusted. However, without permission of the authors, it is not allowed to publish or distribute these files. 

This research was supported by the European Research Council (ERC) under the European Unions Horizon 2020 research and innovation programme (TRUST-949796).



%%%=======================================================================================

-References


[1] M. R. Sheriff, F. F. Redel, and P. Mohajerin Esfahani, Fast Algorithm for Constrained Linear Inverse Problems. arXiv: 2212.01068.

[2] M. V. Afonso, J. M. Bioucas-Dias, and M. A. T. Figueiredo, An Augmented Lagrangian Approach to the Constrained Optimization Formulation of Imaging Inverse Problems, IEEE Transactions On Image Processing, (2009).
    
[3] A. Chambolle and T. Pock. On the Ergodic Convergence rates of a First-order Primal–Dual Algorithm. Mathematical Programming, 159(1-2):253–287, 9 2016
    
[4] J. Duchi, S.S. Shwartz, and T. Chandra. Efficient Projections onto the l1-Ball for Learning in High Dimensions, Google, Technical report, Proceedings of the 25th International Conference on Machine Learning, 2008, Mountain View, CA 94043 




