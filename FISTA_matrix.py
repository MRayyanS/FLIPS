import time
from math import sqrt
import numpy as np
from scipy import linalg


# Soft-thresholding function
def soft_threshold(Z, theta):
    return np.sign(Z) * np.maximum(np.abs(Z) - theta, 0.)

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
def FISTA(X, phi, phiadj_X, phiadj_phi, reg_parameter, GD_stepsize, maxiter, opt_threshold):
    F = np.zeros_like(phiadj_X)
    
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

        obj = np.linalg.norm(X - phi@F, 'fro') ** 2 + reg_parameter * np.sum( np.abs(F) )
        # print(f"iter = {iter}, obj_func = {obj}, opt_check = {opt_check}")

    return F, conv_iter


# generate problem data

# scaling the measurements to satisfy SNR value
SNR         = 10  # in absolute value
sigma_noise = 1   # variance of each entry of noise

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
phiadjX    = phi.T @ X_noisy

algo = 'FISTA'

L = np.max(np.linalg.eig(phiadjphi)[0]) 

GD_stepsize   = 1/L
maxiter       = 500
opt_threshold = 0.0001
F_ISTA        = np.zeros_like(F_true)

# run ISTA with the current reg_parameter
reg_parameter  = 15

if algo == 'ISTA':
    F_rec, conv_iter = ISTA(X_noisy, phi, phiadjX, phiadjphi, reg_parameter, GD_stepsize, maxiter, opt_threshold)
elif algo == 'FISTA':
    F_rec, conv_iter = FISTA(X_noisy, phi, phiadjX, phiadjphi, reg_parameter, GD_stepsize, maxiter, opt_threshold)

# recover X_true from X_noisy 
X_rec = phi @ F_rec

# define error
numer = np.linalg.norm(F_true - F_rec, 'fro') ** 2
denom = np.linalg.norm(F_true, 'fro') ** 2
NMSE  = numer / denom

print(f"|noise| = {np.linalg.norm(noise, 'fro')}, |X - Xrec| = {np.linalg.norm(X_noisy - X_rec, 'fro')}, convereged at iter = {conv_iter}, NMSE = {NMSE}")
