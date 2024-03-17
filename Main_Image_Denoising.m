%% For Denoising
%  These files are alowed to be adjusted. However, without permission of
%  the authors, it is not allowed to publish or distrubute these files.

% clc
clear all
close all 

%% Initialize and import 

% Import image
im.original    = imread('ASMLhq.png');
patchsize      = 8;                                                             
rescale_min    = 0;                                                   
rescale_max    = 1;                                                        
im.original    = rescale(im.original,rescale_min,rescale_max);           

% Selecting image dimensions
bigimdim1    = 930;                                                     
bigimdim2    = 2642;     
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
indices_store_patches_inp= reshape(1:imdim1*imdim2,[imdim1 imdim2]);        

% Extract patches
X.noise      = double(im2col(im.noise,   [patchsize patchsize],'Sliding'));            
X.no_noise   = double(im2col(im.original,[patchsize patchsize],'Sliding'));              
store_inp    = double(im2col(indices_store_patches_inp, [patchsize patchsize],'Sliding'));                         
                                           

% Dimensions 

% n = dimension of each patch = patchsize^2     
% N = no. of patches
[n, N]  = size(X.no_noise);    

% select patch dimensions
patch_dimensions  =  [patchsize patchsize] ;
                                                                                                                                                     

% Selecting Dictionary
dictionary        =  'dct' ;


% Calculating appriopriate epsilon
epsilon = 1*sqrt(var_noise)*sqrt(n);  

maxiter        = 100 ;


% Select appropriate oracle and oracle-parameters

% oracle         = 'SimpleQO' ;
oracle         = 'AcceleratedQO' ;
                                                             
momentum_para  = -10 ;

% Select appropriate smoothness parameter for current patch size 
% (tuned values for cameraman image)
if patchsize ==4
    betainv = 0.975;  
elseif patchsize ==8
    betainv = 0.9;  
elseif patchsize == 16
    betainv = 0.75; 
elseif patchsize ==32
    betainv = 0.475; 
elseif patchsize == 64
    betainv = 0.03; 
elseif patchsize == 128
    betainv = 0.0015; 
elseif patchsize == 256
    betainv = 0.0002 ;
elseif patchsize == 512
    betainv = 0.000025 ;
end


% select first-order optimality threshold = 0.5*( 1 + (eta_grad'*h)/(norm(eta_grad,inf)) ) ;
% it must be a value in [0, 1], and 0 at optimal solution

opt_threshold   = 0.01 ; 


%% FLIPS solver fucntion
betainv
t_FLIPS  = cputime ;
[F, e_flag, iter_till_conv] = FLIPS_denoising_final(X.noise,dictionary,epsilon,maxiter,oracle,betainv,momentum_para,patch_dimensions,opt_threshold) ;
t_FLIPS  = cputime - t_FLIPS


% run 'show_results.m' to recover images and to show images and plot relevant things



