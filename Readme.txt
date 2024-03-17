Read.me FLIPS

'Main' file contains the general example/problem of this package.
Min_f       c(f)
s.t.        ||x – phi(f)||_2 <= epsilon

The problem is solved with the oracles and state-of-the-art algorithms presented in the paper.
The other files contain the following information:
	
- C_SALSA.m --> Contains the C-SALSA solver as described in:
	Manya V. Afonso, José M. Bioucas-Dias, and Mário A. T. Figueiredo. An Augmented
	Lagrangian Approach to the Constrained Optimization Formulation of Imaging Inverse
	Problems. IEEE Transactions On Image Processing, 2009
	


- ChambollePock.m --> Contains the Chambolle-Pock solver as described in:
	Antonin Chambolle and Thomas Pock. On the ergodic convergence rates of a first-order
	primal–dual algorithm. Mathematical Programming, 159(1-2):253–287, 9 2016
	
	
	
- FLIPS_Solver --> Solver FLIPS for the Quadratic Oracle
	
	

These files are allowed to be adjusted. However, without permission of the authors, it is not allowed to publish or distribute these files. 

This research was supported by the European Research Council (ERC) under the European Unions Horizon 2020 research and innovation programme (TRUST-949796).



