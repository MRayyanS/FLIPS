FLIPS MATLAB Package - By M. R. Sheriff, F. F. redel, and P. Mohajerin Esfahani

FLIPS is an algorithm to solve the constrained linear inverse problems (LIP), see [1] for a detailed description of FLIPSS. In particular, the software here is a MATLAB package for the image denoising problem that solves

min_f       ||f||_1 
s.t.        ||x – Df||_2 <= epsilon


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


[1] M. R. Sheriff, F. F. redel, and P. Mohajerin Esfahani, Fast Algorithm for Constrained Linear Inverse Problems. arXiv: 2212.01068.

[2] M. V. Afonso, J. M. Bioucas-Dias, and M. A. T. Figueiredo, An Augmented Lagrangian Approach to the Constrained
    Optimization Formulation of Imaging Inverse Problems, IEEE Transactions On Image Processing, (2009).
    
[3] A. Chambolle and T. Pock. On the ergodic convergence rates of a first-order
	  primal–dual algorithm. Mathematical Programming, 159(1-2):253–287, 9 2016
    
[4] J. Duchi, S.S. Shwartz, and T. Chandra. Efficient Projections onto the l1-Ball for Learning in High Dimensions Google, Technical report, Proceedings of the 25th International Conference on Machine Learning, 2008, Mountain View, CA 94043 




