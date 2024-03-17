function [F, e_flag, iter_till_conv] = FLIPS_denoising_final(X,dictionary,epsilon,maxiter,oracle,betainv,momentum_para,patch_dimensions,opt_threshold)

%% getting dimensions

[n, N] = size(X) ;

%% selecting normalization constant: tau, for well-conditioning
tau = 10 ;
% tau = (1/(N*sqrt(K)))*norm(X,2) ;

% we shall use this to change the constraint on h 
% as ||h||_c \leq tau, instead of ||h||_c \leq 1


iter_till_conv = maxiter*ones(1,N) ;

%% Initialization and storing variables

% computing phi_adj_X

phiadjX    = phi_adjoint_X(X, dictionary, patch_dimensions) ;

% Initializing output quantities

F       = zeros(n,N) ;
e_flag  = zeros(1,N) ;


%% FLIPS-Solver start

for l = 1:N

    % Initialiization
    phiadj_x         = phiadjX(:,l) ;
    normx_minus_eps  = norm(phiadj_x,2)^2 - epsilon^2 ;

% checking if normx is smaller than epsilon assumption
if normx_minus_eps <= 0

        F(:,l)    = zeros(n,1) ;
        e_flag(l) = opt_threshold ;

else

    % initial iterate = l1-normalized least squares solution
    h = (tau/norm(phiadj_x,1))*phiadj_x;

    norm_phih_sq        = h'*h ;
    ip_xphih            = phiadj_x'*h ;

    d_old = zeros(n,1) ;  % only used in accelerated quadratic oracle


    %% iterations start
    for t = 1:maxiter

        %% computing eta function

        % computing the term in the square root
        sqrt_term = ip_xphih^2 - ( norm_phih_sq*normx_minus_eps ) ;
        sqrt_term = sqrt(sqrt_term) ;

        numer = normx_minus_eps ;
        denom = ip_xphih + sqrt_term ;
        eta_val = numer/denom ;

        %% computing grad eta
        
        eta_grad   = eta_val*h - phiadj_x ;
        alpha        = eta_val / sqrt_term ;
        eta_grad   = alpha*eta_grad ;

        % grad_check1 = -eta_grad'*phiadj_x 

        %% Descent direction oracles

        if strcmp(oracle, 'SimpleQO')
            % simple quadratic oracle
            g_oracle_out = Simple_quad_oracle(h, eta_grad, betainv,tau) ;
        elseif strcmp(oracle, 'AcceleratedQO')
            % accelerated quadratic oracle 
            [g_oracle_out, d_old] = Accelerated_quad_oracle(h, eta_grad, d_old, betainv, momentum_para,tau) ;
        end

% phiadjphi_g  = g_oracle_out ; % but not needed, can directly use g_oracle_out
        
%% Exact line search oracle of FLIPS

        ip_xphig     = phiadj_x'*g_oracle_out ;
        ip_xphid     = ip_xphig  -  ip_xphih ;

        norm_phig_sq = g_oracle_out'*g_oracle_out ; % g_oracle_out is phiadjphi_g

        ip_phih_phig = g_oracle_out'*h ;
        ip_phih_phid = ip_phih_phig - norm_phih_sq ;
        ip_phig_phid = norm_phig_sq - ip_phih_phig ;
 
        norm_phid_sq = norm_phig_sq + norm_phih_sq - 2*ip_phih_phig ;


        % computing exact line search

        % quantity to check for gamma = 0
        gamma0_check = ip_xphid - eta_val*ip_phih_phid ;

        % quantities to check for gamma = 1
        % needed to compute eta(g) to see if gamma = 1 is valid
        sqrt_term_g = ip_xphig^2 - normx_minus_eps*norm_phig_sq ;

        % computing eta(g)
        eta_g = normx_minus_eps / ( ip_xphig + sqrt(sqrt_term_g) ) ;

        % condition for gamma = 1
        gamma1_check = ip_xphid - eta_g*ip_phig_phid ;

        if gamma0_check <= 0
            step_size = 0 ;

        % condition for g in cone and step-size = 1
        elseif (sqrt_term_g >= 0) && ( gamma1_check >= 0 )
            step_size = 1 ;
        
        % then gamma \in (0,1) and obtained by solving the quadratic equation
        else 
            % step-size is the positive root of the quadratic eqn

            %setting up the scalar quadratic equation parameters
            a = normx_minus_eps*norm_phid_sq - ip_xphid^2 ;
            b = 2*( normx_minus_eps*ip_phih_phid  -  ip_xphih*ip_xphid ) ;
       
            term1 = normx_minus_eps*ip_phih_phid^2 ;
            term2 = 2*ip_xphid*ip_xphih*ip_phih_phid ;
            term3 = norm_phih_sq*ip_xphid^2 ;

            c = term1 - term2 + term3 ;
            c = c/norm_phid_sq ;

            
            % finding the correct root of the quadratic equation
            root1 = -b + sqrt(b^2 - 4*a*c) ;
            root1 = root1/(2*a) ;

            ip_phid_phihgamma_1 = ip_phih_phid + root1*norm_phid_sq ;
            ip_xphihgamma_1     = ip_xphih + root1*ip_xphid ;

            root1_check = (normx_minus_eps*ip_phid_phihgamma_1/ip_xphid) - ip_xphihgamma_1 ;
        
            root2 = -b - sqrt(b^2 - 4*a*c) ;
            root2 = root2/(2*a) ;

            % to select the correct root
            if root1*root2 >= 0
                if (root1_check >= 0)
                    step_size = root1 ;
                else
                    step_size = root2 ;
                end
            else
                step_size = max(root1 , root2) ;
            end
   
        end

%% Update

        h                 = h + step_size*( g_oracle_out - h ) ;

        ip_xphih          = ip_xphih + step_size*(ip_xphig - ip_xphih) ;

        norm_phih_sq      = norm_phih_sq + step_size^2*norm_phid_sq + 2*step_size*ip_phih_phid ;


        % checking first-order optimality for convergence
        opt_check   =  0.5*( 1 + (eta_grad'*h)/(tau*norm(eta_grad,inf)) ) ;

        if (opt_check <= opt_threshold)
            iter_till_conv(l) = t ;
            break % this stops the algorithm for this patch
        end

    % iteration over t ends
    end


    % define the output quantities for this patch
    % F(:,l)     = ( (phiadj_x'*h)/(norm(h,2))^2 ).*h ; % messes up with
    % conv_threshold criteria for CPU times

    F(:,l)     = eta_val*h ;
    e_flag(l)  = opt_check ;
    

% ending the if statement
end

% ending for loop for l, i.e., for current patch
end

% ending the FLIP-Solver function
end


%% Different descent direction Oracles in FLIPS

% Simple Quadratic Oracle
function [g] = Simple_quad_oracle(h, grad, betainv,tau)

% gradient step for any LIP
g = h - betainv*grad ;

% LIP specific oracle, for sparse regression problems where c = l1-norm

g = Projection_l1_ball(g, tau); % computes projection onto l1-balls of radius tau

end


% Accelerated Quadratic Oracle
function [g, d_new] = Accelerated_quad_oracle(h, grad, d_old, betainv, momentum_para,tau) 

% gradient step for any LIP
g = h - betainv*grad - (betainv*momentum_para)*d_old ;

% LIP specific oracle, for sparse regression problems where c = l1-norm
g = Projection_l1_ball(g, tau) ; % computes projection onto l1-balls of radius tau
d_new = g - h ;

end


% Projection onto l1-ball

function w = Projection_l1_ball(v, tau)
% w = ProjectOntoL1Ball(v, b) returns the vector w which is the solution
%   to the following constrained minimization problem:
%
%    min   ||w - v||_2
%    s.t.  ||w||_1 <= tau.
%
%   It performs Euclidean projection of v to the l1-ball of radius tau.

n = size(v) ;

if (tau < 0)
  error('Radius of L1 ball is negative: %2.3f\n', tau);
elseif (norm(v, 1) <= tau)
  norm(v, 1)
  w = v;
  return;
end

u = sort(abs(v),'descend') ;
u = [u; 0] ; % just helps when the ||v||_1 is already close to tau

cumsum = 0 ;

for i = 1:n
    cumsum = cumsum + u(i) ;
    if u(i+1) < u(i)
        % check if the thresholding value is met
        if cumsum >= tau + i*u(i+1)
            theta = (1/i)*(cumsum - tau) ;

            % do soft thresholding with theta
            w = sign(v).* max(abs(v) - theta, 0) ;
            return;
        end
    end
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


