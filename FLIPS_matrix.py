import time
from math import sqrt
import numpy as np
from scipy import linalg

from utils import *


# generate problem data

# scaling the measurements to satisfy SNR value
SNR         = 20  # in absolute value
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
phiadjX    = phi.T @ X_noisy

# initialization of iterates
tau      =  10000      # constant only relevant to make the problem well-conditioned

F_l2  = np.linalg.lstsq(phi, X_noisy, rcond=None)[0] 
H     = (tau / np.sum(np.abs(F_l2))) * F_l2

# define aparameters of the algorithm
maxiter       = 1000
oracle        = 'SimpleQO'
GD_stepsize   = 1000
momentum_para = 0.9
opt_threshold = 0.00000075


# for current epsilon, run FLIPS
epsilon         = 0.85 * sigma_w * np.sqrt(m * N)
normX_minus_eps = np.linalg.norm(X_noisy, 'fro')**2 - epsilon**2

if normX_minus_eps > 0:
    F_FLIPS, H, conv_FLIPS, avg_stepsize = FLIPS_Solver(phiadjX, phiadjphi, normX_minus_eps, H, maxiter, oracle, GD_stepsize, momentum_para, opt_threshold)
else:
    F_FLIPS    = np.zeros_like(H)
    H          = np.zeros_like(H)
    conv_FLIPS = 0

# recover X_true from X_noisy 
X_FLIPS = phi @ F_FLIPS

# define error
numer = np.linalg.norm(F_true - F_FLIPS, 'fro') ** 2
denom = np.linalg.norm(F_true, 'fro') ** 2
NMSE_FLIPS = numer / denom

print(f"|X - Xrec| = {np.linalg.norm(X_noisy - X_FLIPS, 'fro')}, converged at iter = {conv_FLIPS}, avg_stepsize = {avg_stepsize}, NMSE_FLIPS = {NMSE_FLIPS}")