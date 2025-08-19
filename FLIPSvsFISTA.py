import time
from math import sqrt
import numpy as np

from scipy import linalg

# import all the required functions for sparse coding
from utils import * 


# generate problem data

m  = 350                        # number of measurements
K  = 500                        # dimension of sparse signal
S  = int(np.ceil(0.15*K))       # sparsity level
N  = 100                        # number of signals to solve the sparse-coding problem

# initialize random number generator
rng  = np.random.RandomState(42)

# define linear measurement matrix phi
phi       = rng.randn(m, K)  # random design
# phi       = np.eye(K)
phiadjphi = phi.T @ phi    # compute phiadjphi

# define smoothness constant
L = np.max(np.linalg.eig(phiadjphi)[0])

# variance of each entry of noise
sigma_noise = 1

# number of SNR values
SNR_min = 5
SNR_max = 20
N_SNR   = 4
SNR_seq = np.linspace(SNR_min, SNR_max, N_SNR)  # in absolute value

# number of grid points for regularization parameter reg_param
reg_param_min = 5
reg_param_max = 25
N_reg_param   = 100
reg_param_seq = np.linspace(reg_param_max, reg_param_min, N_reg_param)

# number of grid points of epsilon
eps_seq_type = 'uniform' # choices = 'uniform' and 'FISTA'

N_eps               = N_reg_param
eps_multiples       = np.linspace(1, 0.5, N_eps)
epsilon_seq_uniform = eps_multiples * ( sigma_noise * np.sqrt(m * N) )


# initialize performance metrics
NMSE_seq_FISTA    = {}
Conv_iter_FISTA   = {}

epsilon_seq = {}

 # initialize performance metrics
Conv_iter_FLIPS   = {}
NMSE_seq_FLIPS    = {}


for SNR in SNR_seq:
    # conduct the experiment for this value of SNR

    # define true F
    F_true = rng.randn(S, N)
    F_true = np.append(F_true, np.zeros((K-S, N)), axis=0)
    F_true = permute_columns(F_true)     # permutes entries in each column

    # F_true scaled to satisfy SNR
    SNR_abs = 10 ** ( SNR / 10 )
    sigma_signal = sigma_noise*np.sqrt(K*m*SNR_abs / (S*np.trace(phiadjphi)))
    F_true       = sigma_signal*F_true

    # obtain clean measurements, i.e., X_true 
    X_true  = np.dot(phi, F_true)

    # generate noise 
    noise  = rng.randn(*X_true.shape)
    noise  = sigma_noise*noise

    # noisy measuremremensts X_noisy
    X_noisy = X_true + noise

    # Computing quantities that are global
    phiadjX    = phi.T @ X_noisy

    # run FISTA at this SNR for various values of reg_param
    # select an algorithm to perform the experiment nd its parameters
    algo          = 'FISTA'
    GD_stepsize   = 1/L
    maxiter       = 500
    opt_threshold = 0.00001
    F             = np.zeros_like(F_true)

    # initialize performance metrics
    NMSE_seq_FISTA[SNR]  = np.array([])
    Conv_iter_FISTA[SNR] = np.array([])
    
    epsilon_seq_FISTA    = np.array([])

    for reg_param in reg_param_seq:
        # run ISTA with the current reg_parameter

        if algo == 'ISTA':
            F, conv_iter = ISTA(X_noisy, phi, phiadjX, phiadjphi, reg_param, GD_stepsize, maxiter, opt_threshold)
        elif algo == 'FISTA':
            F, conv_iter = FISTA(phiadjX, phiadjphi, F, reg_param, GD_stepsize, maxiter, opt_threshold)

        # recovered signal from FISTA
        X_rec = phi @ F

        # compute F_FISTA adjusting to scalar scaling
        numer = np.trace( phiadjX.T @ F )
        denom = np.linalg.norm( X_rec, 'fro' ) ** 2
        F_FISTA = (numer/denom) * F

        # compute NMSE
        numer        = np.linalg.norm(F_true - F_FISTA, 'fro') ** 2
        denom        = np.linalg.norm(F_true, 'fro') ** 2
        current_NMSE = 10*np.log10(numer / denom)
        NMSE_seq_FISTA[SNR]  = np.append(NMSE_seq_FISTA[SNR], current_NMSE ) 

        # collect iterations required for convergence
        Conv_iter_FISTA[SNR] = np.append( Conv_iter_FISTA[SNR], conv_iter )

        # compuet the value of epsilon for this experiment
        current_eps = np.linalg.norm(X_noisy - X_rec, 'fro')
        epsilon_seq_FISTA = np.append(epsilon_seq_FISTA, current_eps )
        
        print(f"SNR = {SNR}dB, algo = FISTA, |X - Xrec| = {current_eps}, NMSE = {current_NMSE}, conv iter = {conv_iter}")


    # run FLIPS for several values of epsilon for this curent SNR

    if eps_seq_type == 'FISTA':
        epsilon_seq[SNR] = np.sort(epsilon_seq_FISTA)[::-1] # to get decreasing sequence of eps
    else:
        epsilon_seq[SNR] = epsilon_seq_uniform

    # initialization of iterates
    tau  =  10000
    # tau : constant only relevant to make the problem well-conditioned
    # larger the value, better the convergence, but we have to tune opt_threshold propoerly

    F_l2 = np.linalg.lstsq(phi, X_noisy, rcond=None)[0] 
    H    = (tau / np.sum(np.abs(F_l2))) * F_l2

    # define aparameters of the algorithm
    maxiter       = 1000
    oracle        = 'SimpleQO'
    GD_stepsize   = 1000 
    momentum_para = 0.9
    opt_threshold = 0.00000075

    # initialize performance metrics
    Conv_iter_FLIPS[SNR] = np.array([])
    NMSE_seq_FLIPS[SNR]  = np.array([])

    # for current SNR level sun FLIPS for several (N_eps many) values of epsilon
    for epsilon in epsilon_seq[SNR]:

        # for current epsilon, run FLIPS
        normX_minus_eps = np.linalg.norm(X_noisy, 'fro')**2 - epsilon**2

        if normX_minus_eps > 0:
            F_FLIPS, H, conv_iter, avg_stepsize = FLIPS_Solver(phiadjX, phiadjphi, normX_minus_eps, H, maxiter, oracle, GD_stepsize, momentum_para, opt_threshold)
        else:
            F_FLIPS   = np.zeros_like(H)
            H         = np.zeros_like(H)
            conv_iter = 0

        # define error
        numer = np.linalg.norm(F_true - F_FLIPS, 'fro') ** 2
        denom = np.linalg.norm(F_true, 'fro') ** 2
        current_NMSE = 10*np.log10(numer / denom)
        NMSE_seq_FLIPS[SNR] = np.append( NMSE_seq_FLIPS[SNR], current_NMSE )

        Conv_iter_FLIPS[SNR] = np.append(Conv_iter_FLIPS[SNR], conv_iter)

        # recover X_true from X_noisy 
        X_FLIPS = phi @ F_FLIPS
        print(f"SNR = {SNR}dB, algo = FLIPS, epsilon = {epsilon}, NMSE = {current_NMSE}, conv iter = {conv_iter}, avg_FW_stepsize = {avg_stepsize}")



# save the simulation data using pickle
import pickle

# Save the variables to a file
with open("FISTA_v_FLIPS_uniform_eps_seq.pkl", "wb") as file:
    pickle.dump((NMSE_seq_FISTA, NMSE_seq_FLIPS, Conv_iter_FISTA, Conv_iter_FLIPS, SNR_seq, reg_param_seq, epsilon_seq, m, N, sigma_noise), file)

print("Variables saved successfully!")


