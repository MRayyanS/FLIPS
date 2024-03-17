%% For Compressed sensing
%  These files are alowed to be adjusted. However, without permission of
%  the authors, it is not allowed to publish or distrubute these files.

clc,
close all 

%% Dimensions of Signals
N  = 1 ;    % number of signals                                          
n  = 5000 ;    % dimension of each signal  

%% Generating the data for the problem

% Getting the true signal of +1, -1 entries

theta_n = n/2 ;
F_true = [ones(1, theta_n) -1*ones(1, n-theta_n)]' ;

% defining measurement matrix (or the linear map phi)

m  = ceil(0.55*n);   % number of measurements             
phi = rand(m, n) - 0.5.*ones(m,n) ;


%% Getting measurements, non-noisy and noisy

% non-noisy measurements
X_no_noise   = phi*F_true ;      

% generating measurement noisev'*ones  
std_dev             = 0.0125;
meas_noise          = std_dev*randn(m,1);   

% getting noisy measurements
X_noise             = X_no_noise + meas_noise ; 


%% Calculating appriopriate epsilon

normx = norm(X_noise,2)

noise_norm = norm(meas_noise, 2)

epsilon = 10*std_dev*sqrt(m) 


%% FLIPS (Quadratic Oracle -- QO)

% Select appropriate oracle and oracle-parameters
% (tuned values for cameraman image)

% oracle         = 'SimpleQO' ;
oracle         = 'AcceleratedQO' ;

betainv        = 0.1 ;                                                              
momentum_para  = -10 ;

maxiter        = 4000 ;

% FLIPS solver
[F, eta, gamma, optimality_check, Plot_signal] = FLIPS_Solver(X_noise,phi,epsilon,maxiter,oracle,betainv,momentum_para) ;


%% Recover results

F_rec  =  F(:,end) ;

%%  =======================================================================


% plot using any format, you can run the file `show_results.m' file
% after this








