function [H,time_count,iter_count] = C_SALSA_CPU_time(X,dictionary,epsilon,maxiter,mu1,patch_dimensions,H_star,conv_threshold)
% C-SALSA

time_count = 0 ;
iter_count = 0 ;

t0 = cputime ;

[n,N] = size(X) ;

phiadjX    = phi_adjoint_X(X, dictionary, patch_dimensions) ;

H  = zeros(size(phiadjX)) ;

tau= 1/mu1;
for c = 1:N
    
    % Initialization for every iteration
    
    % working with phit x instead of x to speed up
    phitx = phiadjX(:,c) ;

    normx_minus_eps = norm(phitx,2)^2 - epsilon^2 ;

    % checking if normx is smaller than epsilon assumption
    if normx_minus_eps <= 0

        H(:,c) = zeros(n,1) ;

    else

    % h     = phitx; 
    v1    = 0*phitx;
    v2    = 0*phitx;
    d1    = 0*phitx;
    d2    = 0*phitx;
    
    psi_function = @(h,tau) soft(h,tau);
  
    h_star = H_star(:, c) ;
    norm_h_star = norm(h_star,2) ;
    
    for t=1:maxiter
        
        r = v1 + d1 + (v2 + d2);
        u = 0.5*r;
        v1= psi_function(u-d1,tau);
        Du = u;
        s = Du-d2;
        if norm(s-phitx,2)> epsilon
            v2= phitx+ epsilon*(s-phitx)/norm(s-phitx,2);
        else
            v2= phitx+ s-phitx;
        end
        d1 = d1 - u  + v1;  % Lagrange variable 1
        d2 = d2 - Du + v2;  % Lagrange variable 2
        

        %% check for stopping criteria, if satisfied, then break out of for loop over t

        t1 = cputime ;
        
        convergence_check = norm(u - h_star,2)/norm_h_star ;
        if (convergence_check <= conv_threshold)
            break % this stops the algorithm for this patch
        end

        time_count = time_count + t1 - t0 ;
        t0 = cputime ;

    % iteration over t ends
    end

    % computes average number of iterations over different patches
    iter_count  =  iter_count + (1/c)*(t - iter_count) ;
    H(:,c) = u;
    
    % ending if statement
    end
    
 % ending the for loop on c, i.e., algorithm for the current patch
 end

end




%% computing phi_adjoint_X

function phiadjX    = phi_adjoint_X(X, dictionary, patch_dimensions) 

[n, N] = size(X) ;

phiadjX  = zeros(n,N) ;

for i = 1:N

    patch        = reshape(X(:,i), patch_dimensions) ; % for fast computation of dct/basis coefficients

    if strcmp(dictionary,'dct')
        dct_coeff    = dct2(patch) ;
        phiadjX(:,i) = reshape(dct_coeff, [n 1]) ;
    end

    

end

end


%%  Function for soft-thresholding


function x_soft = soft(x_orig,tau)
% Soft thresholding algorithm 

if sum(abs(tau(:)))~=0
    threshold = abs(x_orig) - tau;
    x_soft = max(threshold, 0);
    x_soft = x_soft./(x_soft + tau) .* x_orig;
else
    x_soft = x_orig;
end

end
