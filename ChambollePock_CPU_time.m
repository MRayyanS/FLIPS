function [H,time_count,iter_count] = ChambollePock_CPU_time(X,dictionary,epsilon,maxiter,tau,sigma,theta,patch_dimensions,H_star,conv_threshold)

% Chambolle-Pock algorithm 

[n,N] = size(X);

time_count = 0 ;
iter_count = 0 ;

t0 = cputime ;

phiadjX    = phi_adjoint_X(X, dictionary, patch_dimensions) ;

H  = zeros(size(phiadjX)) ;

for c = 1:N

    % phiadj x instead of as it speeds up
    phitx     = phiadjX(:,c) ;

    normx_minus_eps     = norm(phitx,2)^2 - epsilon^2 ;

    % checking if normx is smaller than epsilon assumption
    if normx_minus_eps <= 0

        H(:,c) = zeros(n,1) ;

    else

    h         = phitx ; 
    h         = (1/norm(h,1))*h ;
    hplus     = 0*h; 
    u         = 0*h;

    h_star = H_star(:, c) ;
    norm_h_star = norm(h_star,2) ;

    for t=1:maxiter
    
        h_bar = h - tau*u;
        for j=1:n
            hplus(j) = soft(h_bar(j),tau);
        end
        h = hplus + theta*(hplus - h);

        u_bar = u + sigma*h - sigma*phitx;
        u = (1 - epsilon*sigma/norm(u_bar,2))*u_bar;


        %% check for stopping criteria, if satisfied, then break out of for loop over t

        t1 = cputime ;
        
        convergence_check = norm(h - h_star,2)/norm_h_star ;
        
        if (convergence_check <= conv_threshold)
            break % this stops the algorithm for this patch
        end

        time_count = time_count + t1 - t0 ;
        t0 = cputime ;

    % iteration over t ends
    end

    % computes average number of iterations over different patches
    iter_count  =  iter_count + (1/c)*(t - iter_count) ;
    H(:,c) = sparse(h);

    % ending if statement
    end

% iteration over c, i.e., patches, ends
end

% ending the CP function
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







