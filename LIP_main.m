%% For Compressed sensing
%  These files are alowed to be adjusted. However, without permission of
%  the authors, it is not allowed to publish or distrubute these files.

clc,
close all 

%% Importing and generating the data for the problem

% Import image
im_original    = imread('lena.png');
% im_original    = rgb2gray(im_original) ;

%Rescaling for appropriate size
patchsize      = 128;                                                             
rescale_min    = 0;                                                    
rescale_max    = 1;                                                        
im_original    = rescale( im_original,rescale_min,rescale_max);    
im_original    = imresize(im_original,[patchsize nan]);
input_im_size  = size(im_original) ;

% Getting a matrix of sequence of patches along its columns

F_true   =  double(im2col(im_original,[patchsize patchsize],'Sliding'));     

%% Dimensions 


% of Signals
N  = size(F_true,2) ;    % number of signals                                          
n  = size(F_true,1) ;    % dimension of each signal  

% of Dictionary
K  = n;                           % No of atoms in the dictionary 
D  = DCT(K);                      % DCT Dictionary

%% defining measurement matrix, and lin-map phi

% measurement vector size

m  = ceil(0.65*n);                          

% Measurement matrix 

C = rand(m, n) - 0.5.*ones(m,n) ;   % for compressed sensing  

% Linear mapping phi in the LIP

phi = C*D; 

% phi = eye(n,n) ; % to do testing

%% Getting measurements, non-noisy and noisy

% non-noisy measurements
X_no_noise   = C*F_true ;      

% generating measurement noise via guassian image
mean_noiselevel     = 0;       
var_noise           = 0.0025;
meas_noise          = imnoise(zeros(size(X_no_noise)),'gaussian',mean_noiselevel,var_noise);

X_noise             = X_no_noise + 0.05*meas_noise ; 

%% Calculating appriopriate epsilon

normx = norm(X_noise,2)

noise_norm = norm(meas_noise, 2)

% epsilon = 1*sqrt(var_noise)*sqrt(m) 

epsilon = 5*noise_norm

%% FLIPS (Quadratic Oracle -- QO)

% Select appropriate oracle and oracle-parametersc
% (tuned values for cameraman image)

% oracle         = 'SimpleQO' ;
oracle         = 'AcceleratedQO' ;

betainv        = 0.000125 ;                                                              
momentum_para  = -500 ;

maxiter        = 750 ;

% FLIPS solver
[F, eta, gamma, Plot_signal] = FLIPS_Solver(X_noise,phi,epsilon,maxiter,oracle,betainv,momentum_para) ;










%% DCT function used in main

function D = DCT(n)
% DCT-dictionary function 

D = zeros(n,n);
for i = 1:n
    e      = zeros(1,n);
    e(i)   = 1;
    D(:,i) = idct2(e);
end

end



