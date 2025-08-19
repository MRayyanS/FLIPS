import time
from math import sqrt
import numpy as np
from scipy import linalg

# Soft-thresholding function
def soft_thresh(f, l):
    return np.sign(f) * np.maximum(np.abs(f) - l, 0.)

## writing FLIPS function
def FLIPS_Solver(phiadjx, phiadjphi, normx_minus_eps, h, maxiter, oracle, GD_stepsize, momentum_para, opt_threshold):
    
    # Selecting normalization constant: tau, for well-conditioning of the problem
    tau = np.linalg.norm(h, 1)

    conv_iter = maxiter

    # initialization of iterating quantities
    phiadjphi_h   = phiadjphi @ h
    norm_phih_sq  = h.T @ phiadjphi_h
    ip_xphih      = h.T @ phiadjx
        
    g_oracle_out  = np.zeros_like(h)
    d_old         = np.zeros_like(h)  # Only used in accelerated quadratic oracle
    
    # FLIPS - iterations start
    for iter in range(maxiter):
        
        # Computing the term in the square root
        sqrt_term = ip_xphih ** 2 - (norm_phih_sq * normx_minus_eps)
        sqrt_term = np.sqrt(sqrt_term)

        # Computing eta function
        numer = normx_minus_eps
        denom = ip_xphih + sqrt_term
        eta_val = numer / denom

        # Computing grad eta
        grad  = eta_val * phiadjphi_h - phiadjx
        alpha = eta_val / sqrt_term
        grad  = alpha * grad

         # checking first-order optimality for convergence
        opt_check = np.linalg.norm(grad, np.inf) + (1/tau) * np.dot(grad, h)
        if (opt_check <= opt_threshold):
            # print('converged at iter = ', iter)
            conv_iter = iter
            break # this stops the algorithm for this sample

        # Descent direction oracles
        if oracle == 'SimpleQO':
            g_oracle_out = Simple_quad_oracle(h, grad, GD_stepsize, tau)
        elif oracle == 'AcceleratedQO':
            g_oracle_out, d_old = Accelerated_quad_oracle(h, grad, d_old, GD_stepsize, momentum_para, tau)
            # the output d_old is actually the direction of current update, but d_old is only used in next iteration whereby it is indeed the direction of previous update

        # costly matrix multiplication
        phiadjphi_g = phiadjphi @ g_oracle_out

        # computing quantitities for exact line search
        ip_xphig = phiadjx.T @ g_oracle_out
        ip_xphid = ip_xphig - ip_xphih

        norm_phig_sq = g_oracle_out.T @ phiadjphi_g
        ip_phih_phig = phiadjphi_h.T @ g_oracle_out
        ip_phih_phid = ip_phih_phig - norm_phih_sq
        ip_phig_phid = norm_phig_sq - ip_phih_phig

        norm_phid_sq = norm_phig_sq + norm_phih_sq - 2 * ip_phih_phig

        # Computing exact line search
        gamma0_check = ip_xphid - ( eta_val * ip_phih_phid )
        sqrt_term_g  = ip_xphig ** 2 - normx_minus_eps * norm_phig_sq

        # checking if G is in the cone and if, then compute condition for step-size = 1
        if sqrt_term_g >= 0:  # G(H) is inside the cone
            eta_g = normx_minus_eps / (ip_xphig + np.sqrt(sqrt_term_g))
            gamma1_check = ip_xphid - (eta_g * ip_phig_phid)
        
        if gamma0_check <= 0:
            step_size = 0
        elif sqrt_term_g >= 0 and gamma1_check >= 0:  # G(H) is inside the cone and step-size = 1
            step_size = 1
        else:
            a = normx_minus_eps * norm_phid_sq - ip_xphid ** 2
            b = 2 * (normx_minus_eps * ip_phih_phid - ip_xphih * ip_xphid)
            term1 = normx_minus_eps * ip_phih_phid ** 2
            term2 = 2 * ip_xphid * ip_xphih * ip_phih_phid
            term3 = norm_phih_sq * ip_xphid ** 2

            c = term1 - term2 + term3
            c = c / norm_phid_sq

            root1 = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
            ip_phid_phihgamma_1 = ip_phih_phid + root1 * norm_phid_sq
            ip_xphihgamma_1 = ip_xphih + root1 * ip_xphid
            root1_check = (normx_minus_eps * ip_phid_phihgamma_1 / ip_xphid) - ip_xphihgamma_1

            root2 = (-b - np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)

            if root1 * root2 >= 0:
                if root1_check >= 0:
                    step_size = root1
                else:
                    step_size = root2
            else:
                step_size = max(root1, root2)

        # FLIPS update
        h            = h            + step_size * (g_oracle_out - h)
        ip_xphih     = ip_xphih     + step_size * (ip_xphig - ip_xphih)
        phiadjphi_h  = phiadjphi_h  + step_size * (phiadjphi_g - phiadjphi_h)
        norm_phih_sq = norm_phih_sq + step_size * 2 * ip_phih_phid + step_size ** 2 * norm_phid_sq
    
    f = eta_val*h
    
    return f, h, conv_iter

## Auxialiary functions needed
def Simple_quad_oracle(h, grad, betainv, tau):
    g = h - betainv * grad
    g = Projection_l1_ball(g, tau)
    return g

def Accelerated_quad_oracle(h, grad, d_old, betainv, momentum_para, tau):
    g = h - betainv * grad - (betainv * momentum_para) * d_old
    g = Projection_l1_ball(g, tau)
    d_new = g - h
    return g, d_new

def Projection_l1_ball(v, tau):
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
                return np.sign(v) * np.maximum(np.abs(v) - theta, 0)

# Function to permute each column individually
def permute_columns(matrix):
    permuted_matrix = np.zeros_like(matrix)
    for col in range(matrix.shape[1]):
        permuted_matrix[:, col] = np.random.permutation(matrix[:, col])
    return permuted_matrix



# generate problem data

# scaling the measurements to satisfy SNR value
SNR         = 10  # in absolute value
sigma_noise = 1   # variance of each entry of noise
# sigma_f is what every entry of F_true is multiplied with according to SNR, and it is the variance of each non-zero entry of F_true

rng  = np.random.RandomState(42)

m  = 300                        # number of measurements
K  = 500                        # dimension of sparse signal
s  = int(np.ceil(0.15*K))       # sparsity level
N  = 100                        # number of signals to solve the sparse-coding problem


# define linear measurement matrix phi
phi       = rng.randn(m, K)  # random design
# phi       = np.eye(K)
phiadjphi = phi.T @ phi    # compute phiadjphi 

# define true F
F_true = rng.randn(s, N)
F_true = np.append(F_true, np.zeros((K-s, N)), axis=0)
F_true = permute_columns(F_true)     # permutes entries in each column

# F_true scaled to satisfy SNR
sigma_signal = sigma_noise*np.sqrt(K*m*SNR / (s*np.trace(phiadjphi)))
F_true       = sigma_signal*F_true

# obtain clean measurements, i.e., X_true 
X_true  = np.dot(phi, F_true)

# generate noise
sigma_w  = 1 
noise    = rng.randn(*X_true.shape)
noise    = sigma_w*noise

# noisy measuremremensts X_noisy
X_noisy = X_true + noise

# Computing quantities that are global
phiadj_X    = phi.T @ X_noisy

# initialization of iterates
tau      =  10      # constant only relevant to make the problem well-conditioned

F_FLIPS     = np.linalg.lstsq(phi, X_noisy, rcond=None)[0] 
H           = np.zeros(F_FLIPS.shape)
conv_FLIPS  = np.zeros((N, 1))
for t in range(N):
    H[:, t] = (tau / np.linalg.norm(F_FLIPS[:, t], 1)) * F_FLIPS[:, t]


# define aparameters of the algorithm
maxiter       = 200
oracle        = 'SimpleQO'
GD_stepsize   = 5
momentum_para = 0.005
opt_threshold = 0.001

for t in range(N):
    x        = X_noisy[:,t] 
    phiadjx  = phiadj_X[:,t]
    h        = H[:,t]
    
    # define epsilon and run FLIPS
    epsilon         = 0.875 * sigma_w * np.sqrt(m)
    normx_minus_eps = np.linalg.norm(x, axis=0)**2 - epsilon**2
    
    if normx_minus_eps > 0:
        F_FLIPS[:,t], H[:,t], conv_FLIPS[t] = FLIPS_Solver(phiadjx, phiadjphi, normx_minus_eps, h, maxiter, oracle, GD_stepsize, momentum_para, opt_threshold)
    else:
        F_FLIPS[:,t]  = np.zeros_like(h)
        H[:,t]        = np.zeros_like(h)
        conv_FLIPS[t] = 0

    print('epsilon = ', epsilon, 'Sample = ', t, '|noise| = ', np.linalg.norm(noise[:,t],2), '|X - Xrec| = ', np.linalg.norm(X_noisy[:,t] - (phi @ F_FLIPS[:,t]), 2), 'converged at iteration = ', conv_FLIPS[t])


# recover X_true from X_noisy 
X_rec = phi @ F_FLIPS

# define error
numer = np.linalg.norm(F_true - F_FLIPS, axis = 0) ** 2
denom = np.linalg.norm(F_true, axis = 0) ** 2
error = numer / denom

# create plots
import matplotlib.pyplot as plt
plt.close('all')

# Plot the error array 
plt.plot(error, label = r'$ \frac{\| F_{true} - F_{opt} \|_2^2}{\| F_{true} \|_2^2} $')

# Add a horizontal line at epsilon
# plt.axhline(y=epsilon, color='r', linestyle='--', label = fr'$\epsilon  = {epsilon} $')

# Add labels and title 
plt.xlabel('Column Index')
plt.title('NMSE values from the true solution')

# Add legend 
plt.legend() 

# Show the plot 
plt.show()




