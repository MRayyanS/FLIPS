import time
from math import sqrt
import numpy as np
from scipy import linalg

## writing FLIPS function
def FLIPS_Solver(phiadjX, phiadjphi, normX_minus_eps, H, maxiter, oracle, GD_stepsize, momentum_para, opt_threshold):
    
    # Selecting normalization constant: tau, for well-conditioning of the problem
    tau = np.sum(np.abs(H))

    # select conv_iter = max_iter
    conv_iter = maxiter
    avg_FW_stepsize = 0

    # initialization of iterating quantities
    phiadjphi_H   = phiadjphi @ H
    norm_phiH_sq  = np.trace( H.T @ phiadjphi_H )
    ip_XphiH      = np.trace( H.T @ phiadjX )
        
    G_oracle_out  = np.zeros_like(H)
    D_old         = np.zeros_like(H)  # Only used in accelerated quadratic oracle
    
    # FLIPS - iterations start
    for iter in range(maxiter):
        
        # Computing the term in the square root
        sqrt_term = ip_XphiH ** 2 - (norm_phiH_sq * normX_minus_eps)
        sqrt_term = np.sqrt(sqrt_term)

        # Computing eta function
        numer = normX_minus_eps
        denom = ip_XphiH + sqrt_term
        eta_val = numer / denom

        # Computing Grad eta
        Grad  = eta_val * phiadjphi_H - phiadjX
        alpha = eta_val / sqrt_term
        Grad  = alpha * Grad

         # checking first-order optimality for convergence
        opt_check = np.max(np.abs(Grad)) + (1/tau) * np.trace(Grad.T @ H)

        if (opt_check <= opt_threshold):
            # print('converged at iter = ', iter)
            conv_iter = iter + 1
            break # this stops the algorithm for this sample

        # Descent direction oracles
        if oracle == 'SimpleQO':
            G_oracle_out = Simple_quad_oracle(H, Grad, GD_stepsize, tau)
        elif oracle == 'AcceleratedQO':
            G_oracle_out, D_old = Accelerated_quad_oracle(H, Grad, D_old, GD_stepsize, momentum_para, tau)
            # the output D_old is actually the direction of current update, but D_old is only used in next iteration whereby it is indeed the direction of previous update

        # costly matrix multiplication
        phiadhphi_G = phiadjphi @ G_oracle_out

        # computing quantitities for exact line search
        ip_XphiG = np.trace( phiadjX.T @ G_oracle_out )
        ip_XphiD = ip_XphiG - ip_XphiH

        norm_phiG_sq = np.trace( G_oracle_out.T @ phiadhphi_G )
        ip_phiH_phiG = np.trace( phiadjphi_H.T @ G_oracle_out )
        ip_phiH_phiD = ip_phiH_phiG - norm_phiH_sq
        ip_phiG_phiD = norm_phiG_sq - ip_phiH_phiG

        norm_phiD_sq = norm_phiG_sq + norm_phiH_sq - 2 * ip_phiH_phiG

        # Computing exact line search
        gamma0_check = ip_XphiD - ( eta_val * ip_phiH_phiD )
        sqrt_term_g  = ip_XphiG ** 2 - normX_minus_eps * norm_phiG_sq

        # checking if G is in the cone and if, then compute condition for step-size = 1
        if sqrt_term_g >= 0:  # G(H) is inside the cone
            eta_G = normX_minus_eps / (ip_XphiG + np.sqrt(sqrt_term_g))
            gamma1_check = ip_XphiD - (eta_G * ip_phiG_phiD)
        
        if gamma0_check <= 0:
            step_size = 0
        elif sqrt_term_g >= 0 and gamma1_check >= 0:  # G(H) is inside the cone and step-size = 1
            step_size = 1
        else:
            a = normX_minus_eps * norm_phiD_sq - ip_XphiD ** 2
            b = 2 * (normX_minus_eps * ip_phiH_phiD - ip_XphiH * ip_XphiD)
            term1 = normX_minus_eps * ip_phiH_phiD ** 2
            term2 = 2 * ip_XphiD * ip_XphiH * ip_phiH_phiD
            term3 = norm_phiH_sq * ip_XphiD ** 2

            c = term1 - term2 + term3
            c = c / norm_phiD_sq

            root1 = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
            ip_phid_phihgamma_1 = ip_phiH_phiD + root1 * norm_phiD_sq
            ip_XphiHgamma_1 = ip_XphiH + root1 * ip_XphiD
            root1_check = (normX_minus_eps * ip_phid_phihgamma_1 / ip_XphiD) - ip_XphiHgamma_1

            root2 = (-b - np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)

            if root1 * root2 >= 0:
                if root1_check >= 0:
                    step_size = root1
                else:
                    step_size = root2
            else:
                step_size = max(root1, root2)

        avg_FW_stepsize = avg_FW_stepsize + step_size

        # FLIPS update
        H            = H            + step_size * (G_oracle_out - H)
        ip_XphiH     = ip_XphiH     + step_size * (ip_XphiG - ip_XphiH)
        phiadjphi_H  = phiadjphi_H  + step_size * (phiadhphi_G - phiadjphi_H)
        norm_phiH_sq = norm_phiH_sq + step_size * 2 * ip_phiH_phiD + step_size ** 2 * norm_phiD_sq

        # print to tune step-size, can be deleted later
        # print('iter = ', iter, 'step-size = ', step_size ,'eta = ', eta_val, 'opt_check = ', opt_check)


        # automatically tune scale of GD_stepsize based on FW_stepsize on first iteration
        if iter == 1:
            GD_stepsize = 100 * GD_stepsize * step_size

    
    # computing output vector F
    numer = ip_XphiH
    denom = norm_phiH_sq
    F     = (numer/denom) * H

    # # remove this if unnormalized
    # F = eta_val * H

    avg_FW_stepsize = avg_FW_stepsize / conv_iter
    
    return F, H, conv_iter, avg_FW_stepsize

## Auxialiary functions needed
def Simple_quad_oracle(H, Grad, GD_stepsize, tau):
    G = H - GD_stepsize * Grad
    G = Projection_l1_ball(G, tau)
    return G

def Accelerated_quad_oracle(H, Grad, D_old, GD_stepsize, momentum_para, tau):
    G = H - GD_stepsize * Grad - (GD_stepsize * momentum_para) * D_old
    G = Projection_l1_ball(G, tau)
    D_new = G - H
    return G, D_new

def Projection_l1_ball(V, tau):
    v = V.flatten()

    if tau < 0:
        raise ValueError(f'Radius of L1 ball is negative: {tau:.3f}')
    elif np.linalg.norm(v, 1) < tau:
        return v

    u = np.sort(np.abs(v))[::-1]
    u = np.append(u, 0)

    cumsum = 0
    for i in range(len(v)):
        cumsum += u[i]
        if u[i + 1] < u[i]:
            if cumsum >= tau + (i + 1) * u[i + 1]:
                theta = (cumsum - tau) / (i + 1)
                return np.sign(V) * np.maximum(np.abs(V) - theta, 0)

# Function to permute each column individually
def permute_columns(matrix):
    permuted_matrix = np.zeros_like(matrix)
    for col in range(matrix.shape[1]):
        permuted_matrix[:, col] = np.random.permutation(matrix[:, col])
    return permuted_matrix


# Soft-thresholding function
def soft_threshold(Z, theta):
    return np.sign(Z) * np.maximum(np.abs(Z) - theta, 0.)


# ISTA function
def ISTA(X, phi, phiadj_X, phiadj_phi, reg_parameter, GD_stepsize, maxiter, opt_threshold):
    F = np.zeros_like(phiadj_X)
    
    # ISTA iterations start
    conv_iter = maxiter
    theta     = GD_stepsize * reg_parameter
    for iter in range(maxiter):

        # Proximal gradient step
        Grad = phiadj_phi @ F - phiadj_X
        G    = soft_threshold( F - GD_stepsize * Grad , theta )

        # check for convergence
        opt_check = np.linalg.norm(F - G, 'fro')
        if opt_check <= opt_threshold:
            conv_iter = iter 
            break

        F = G     # G is used as a separate variable to compute convergence criteria

        # print('iter = ', iter, 'obj_func =',  np.linalg.norm(X - phi @ F, 'fro' ) ** 2 + reg_parameter * np.sum(np.abs(F)), 'opt_check = ', opt_check, '|Grad| = ', np.linalg.norm(Grad, 'fro'))


    return F, conv_iter

# FISTA function
def FISTA(phiadj_X, phiadj_phi, F0, reg_parameter, GD_stepsize, maxiter, opt_threshold):
    
    F = F0 # initialize F
    
    # FISTA iterations start
    conv_iter = maxiter
    theta     = GD_stepsize * reg_parameter
    t_old     = 1
    G_old     = F
    for iter in range(maxiter):
        
        # Proximal gradient step
        Grad = phiadj_phi @ F - phiadj_X
        G    = soft_threshold( F - GD_stepsize * Grad , theta )
        
        # check for convergence
        opt_check = np.linalg.norm(F - G, 'fro')
        if opt_check <= opt_threshold:
            conv_iter = iter 
            break

        t = 0.5 + np.sqrt(0.25 + t_old ** 2)
        
        # update
        F     = G + ( (t_old - 1)/t ) * (G - G_old)
        G_old = G
        t_old = t

    return F, conv_iter


