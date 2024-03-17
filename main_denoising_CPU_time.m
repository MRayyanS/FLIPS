%% FLIPS
%  These files are alowed to be adjusted. However, without permission of
%  the authors, it is not allowed to publish or distrubute these files.

% clc, 
% clear all
% close all 

%% Initialize and import 

% Import image
im.original    = imread('cameraman.png');
patchsize      = 128;                                                             
rescale_min    = 0;                                                    
rescale_max    = 1;                                                        
im.original    = rescale( im.original,rescale_min,rescale_max);           

% Selecting image dimensions
bigimdim1    = 128;                                                     
bigimdim2    = 128;     
var_noise    = 0.0055;
start        = 1;                       % Starting pixel
endd         = start + bigimdim1 - 1;   % Final pixel 

% Resize image and adding noise 
im.original              = imresize(im.original,[bigimdim1 bigimdim2]);                   
im.original              = im.original(start:endd,start:endd);
[imdim1, imdim2]         = size(im.original);
rng(1)                                                             
mean_noiselevel          = 0;                                                    
im.noise                 = imnoise(im.original,'gaussian',mean_noiselevel,var_noise);        

% Extract patches
X.noise      = double(im2col(im.noise,   [patchsize patchsize],'Sliding'));            
X.no_noise   = double(im2col(im.original,[patchsize patchsize],'Sliding')); 


indices_store_patches_inp  = reshape(1:imdim1*imdim2,[imdim1 imdim2]); 
store_inp                  = double(im2col(indices_store_patches_inp, [patchsize patchsize],'Sliding'));    % needed for recovery

                                           
% Dimensions 
N  = size(X.no_noise,2);                                               
n  = size(X.no_noise,1);                                            
k  = n;                                                                                                                                                     

%% Calculating appriopriate epsilon and other global parameters

epsilon = sqrt(var_noise)*sqrt(n);   

dictionary  = 'dct' ;

patch_dimensions  =  [patchsize patchsize] ;


%% Calculating optimal solution using FLIPS with large number of iteration

maxiter_conv  = 100 ;

% oracle         = 'SimpleQO' ;
oracle         = 'AcceleratedQO' ;
                                                           
momentum_para  = -10 ;

% Select appropriate smoothness parameter for current patch size 
% (tuned values for cameraman image)
if patchsize == 4
    betainv = 2.5;  
elseif patchsize == 8
    betainv = 1.75;  
elseif patchsize == 16
    betainv = 1.25; 
elseif patchsize == 32
    betainv = 0.4; % 0.45 for lena and barbara, 
elseif patchsize == 64
    betainv = 0.075; 
elseif patchsize == 128
    betainv = 0.0015; 
elseif patchsize == 256
    betainv = 0.0002 ;
elseif patchsize == 512
    betainv = 0.000025 ;
end


% select first-order optimality threshold = 0.5*( 1 + (eta_grad'*h)/(norm(eta_grad,inf)) ) ;
% it must be a value in [0, 1], and 0 at optimal solution

opt_threshold   = 0.001 ;


F_star = FLIPS_denoising_final(X.noise,dictionary,epsilon,maxiter_conv,oracle,betainv,momentum_para,patch_dimensions,opt_threshold) ;


%% setting threshold for covergence of FIPS, CP, and C-SALSA in computing CPU times

conv_threshold = 0.01  

%% FLIPS (Quadratic Oracle)

'FLIPS'

maxiter_FLIPS  = 100 ;

[F_FLIPS, time_FLIPS, iter_FLIPS] = FLIPS_CPU_time(X.noise,dictionary,epsilon,maxiter_FLIPS,oracle,betainv,momentum_para,patch_dimensions,F_star,conv_threshold) ;

betainv
time_FLIPS
iter_FLIPS

%% Chambolle-Pock

 'CP'

L_CP       = 1 ;
tau        = 1/L_CP;
sigma      = 1/L_CP;
theta      = 0.8;

maxiter_CP = 50;

[F_CP, time_CP, iter_CP] = ChambollePock_CPU_time(X.noise,dictionary,epsilon,maxiter_CP,tau,sigma,theta,patch_dimensions,F_star,conv_threshold);

time_CP
iter_CP


%% C-SALSA 

'CSALSA'

maxiter_csalsa = 50;

% Selecting tuned values for mu (tuned for cameraman image)
if     patchsize == 4 
       mu = 2.5 ;
elseif patchsize == 8 
       mu = 5;    
elseif patchsize == 16
       mu = 5 ;
elseif patchsize == 32
       mu = 5 ;
elseif patchsize == 64 
       mu = 8 ;
elseif patchsize == 128
       mu = 8;  
elseif patchsize == 256 
       mu = 9 ;
elseif patchsize == 512
       mu = 15;
end

[F_csalsa, time_csalsa, iter_csalsa] = C_SALSA_CPU_time(X.noise,dictionary,epsilon,maxiter_csalsa,mu,patch_dimensions,F_star,conv_threshold) ;

time_csalsa
iter_csalsa

% Recovering output images by running show_results.m






