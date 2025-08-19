import time
from math import sqrt
import numpy as np
from scipy import linalg




# Soft-thresholding function
def soft_threshold(f, reg_parameter):
    return np.sign(f) * np.maximum(np.abs(f) - reg_parameter, 0.)


# Function to permute each column individually
def permute_columns(matrix):
    permuted_matrix = np.zeros_like(matrix)
    for col in range(matrix.shape[1]):
        permuted_matrix[:, col] = np.random.permutation(matrix[:, col])
    return permuted_matrix



# ISTA function
def ista(phi, x, reg_parameter, maxiter, opt_threshold):
    f = np.zeros(phi.shape[1])
    pobj = []
    L = linalg.norm(phi) ** 2  # Lipschitz constant
    time0 = time.time()

    conv_iter = maxiter
    for iter in range(maxiter):
        g = soft_threshold(f + np.dot(phi.T, x - phi.dot(f)) / L, reg_parameter / L)

        opt_check = np.linalg.norm(f - g, 2)
        
        # check for convergence
        if opt_check <= opt_threshold:
            conv_iter = iter 
            break

        f = g
        this_pobj = 0.5 * linalg.norm(phi.dot(f) - x) ** 2 + reg_parameter * linalg.norm(f, 1)
        pobj.append((time.time() - time0, this_pobj))

    times, pobj = map(np.array, zip(*pobj))
    return f, pobj, times, conv_iter



# FISTA function
def fista(phi, x, reg_parameter, maxit):
    f = np.zeros(phi.shape[1])
    pobj = []
    t = 1
    z = f.copy()
    L = linalg.norm(phi) ** 2
    time0 = time.time()
    for _ in range(maxit):
        f_old = f.copy()
        z = z + phi.T.dot(x - phi.dot(z)) / L
        f = soft_threshold(z, reg_parameter / L)
        t0 = t
        t = (1. + sqrt(1. + 4. * t ** 2)) / 2.
        z = f + ((t0 - 1.) / t) * (f - f_old)
        this_pobj = 0.5 * linalg.norm(phi.dot(f) - x) ** 2 + reg_parameter * linalg.norm(f, 1)
        pobj.append((time.time() - time0, this_pobj))

    times, pobj = map(np.array, zip(*pobj))
    return f, pobj, times




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




maxiter = 1000
opt_threshold = 0.001
F_ISTA    = np.zeros_like(F_true)
conv_ISTA = np.zeros((N,1))

for t in range(N):
    x = X_noisy[:,t]
    reg_parameter  = 12
    
    F_ISTA[:,t], pobj_ISTA, times_ISTA, conv_ISTA[t] = ista(phi, x, reg_parameter, maxiter, opt_threshold)

    print('Sample = ', t, '|noise| = ', np.linalg.norm(noise[:,t],2), '|X - Xrec| = ', np.linalg.norm(X_noisy[:,t] - (phi @ F_ISTA[:,t]), 2), 'convereged at iter = ', conv_ISTA[t])

# recover X_true from X_noisy 
X_ISTA = phi @ F_ISTA

# define error
numer = np.linalg.norm(F_true - F_ISTA, axis = 0) ** 2
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