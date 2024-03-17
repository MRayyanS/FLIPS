function [F, eta, gamma, optimality_check, Plot_signal] = FLIPS_Solver(X,phi,epsilon,maxiter,oracle,betainv,momentum_para)

%% getting dimensions

K = size(phi,2) ;
N = size(X,2) ;

% randomly selecting the signal for plotting purpose
Plot_signal = 1 ;

%% selecting normalization constant: tau, for well-conditioning
tau = 1 ;
% tau = (1/(N*sqrt(K)))*norm(X,2) ;

% we shall use this to change the constraint on h 
% as ||h||_c \leq tau, instead of ||h||_c \leq 1

%% Initial iterate = l1-normalized least squares solution

F0 = phi\X;
H  = zeros(size(F0)) ;

for l=1:N
    H(:,l) = (tau/norm(F0(:,l),1))*F0(:,l);
end

%% Initialization for storing variables

eta         = zeros(1,maxiter) ;
gamma       = zeros(1,maxiter) ;
F           = zeros(K,maxiter) ;

%% First-order optimality quantification

optimality_check = zeros(1, maxiter) ;

%% global constants

% Computing the required quantities

normx_minus_eps    = zeros(1,N) ;
norm_phiH_sq       = zeros(1,N) ;
ip_xphiH           = zeros(1,N) ;

phiadjX            = phi'*X ;
phiadjphi_H        = phi'*phi*H ;

for l = 1:N
    h = H(:,l) ;
    x = X(:,l) ;

    normx_minus_eps(l)     = norm(x,2)^2 - epsilon^2 ;
    norm_phiH_sq(l)        = h'*phiadjphi_H(:,l) ;
    ip_xphiH(l)            = phiadjX(:,l)'*h ;
end

eta_val      = zeros(1,N) ;
eta_grad     = zeros(K,N) ;

step_size    = zeros(1,N) ;

G_oracle_out = zeros(K,N) ;
D_old        = zeros(K,N) ; % only used in accelerated quadratic oracle


%% FLIPS - iterations start

for t = 1:maxiter
    
    for l = 1:N

        % computing the term in the square root
        sqrt_term = ip_xphiH(l)^2 - ( norm_phiH_sq(l)*normx_minus_eps(l) ) ;
        sqrt_term = sqrt(sqrt_term) ;

        
        % computing eta function

        numer = normx_minus_eps(l) ;
        denom = ip_xphiH(l) + sqrt_term ;
        eta_dummy = numer/denom ;

        eta_val(l) = eta_dummy ;


        % computing grad eta
        
        grad_dummy = eta_val(l)*phiadjphi_H(:,l) - phiadjX(:,l) ;
        alpha = eta_dummy / sqrt_term ;
        grad_dummy = alpha*grad_dummy ;

        eta_grad(:,l) = grad_dummy ;


        %% Descent direction oracles

        if strcmp(oracle, 'SimpleQO')
            % simple quadratic oracle
            g = Simple_quad_oracle(H(:,l), eta_grad(:,l), betainv,tau) ;
        elseif strcmp(oracle, 'AcceleratedQO')
            % accelerated quadratic oracle 
            d_old = D_old(:,l) ;
            [g, D_old(:,l)] = Accelerated_quad_oracle(H(:,l), eta_grad(:,l), d_old, betainv, momentum_para,tau) ;
        end

        G_oracle_out(:,l) = g ;
        
    end

    % Multiplying the grammian operator and the oracle output G,
    % computationally efficient if done like this

    % most time consuming operation of the algorithm, can be done
    % alternatively if easier

    phiadjphi_G = phi'*phi*G_oracle_out ;


%% Exact line search and update

    for l = 1:N

%% Exact line search oracle of FLIPS


        ip_xphig     = phiadjX(:,l)'*G_oracle_out(:,l) ;
        ip_xphid     = ip_xphig  -  ip_xphiH(l) ;

        norm_phig_sq = G_oracle_out(:,l)'*phiadjphi_G(:,l) ;

        ip_phih_phig = G_oracle_out(:,l)'*phiadjphi_H(:,l) ;
        ip_phih_phid = ip_phih_phig - norm_phiH_sq(l) ;
        ip_phig_phid = norm_phig_sq - ip_phih_phig ;
 
        norm_phid_sq = norm_phig_sq + norm_phiH_sq(l) - 2*ip_phih_phig ;

% computing optimal step-size

        % quantity to check for gamma = 0
        gamma0_check = ip_xphid - eta_val(l)*ip_phih_phid ;


        % quantities to check for gamma = 1

        % needed to compute eta(g) to see if gamma = 1 is valid
        sqrt_term_g = ip_xphig^2 - normx_minus_eps(l)*norm_phig_sq ;

        % computing eta(g)
        eta_g = normx_minus_eps(l) / ( ip_xphig + sqrt(sqrt_term_g) ) ;

        % condition for gamma = 1
        gamma1_check = ip_xphid - eta_g*ip_phig_phid ;

        if gamma0_check <= 0
            step_size(l) = 0 ;

        % condition for g in cone and step-size = 1
        elseif (sqrt_term_g >= 0) && ( gamma1_check >= 0 )
            step_size(l) = 1 ;
        
        % then gamma \in (0,1) and obtained by solving the quadratic equation
        else 
            % step-size is the root of the quadratic eqn that satisfies
            % root_check condition

            %setting up the scalar quadratic equation parameters
            a = normx_minus_eps(l)*norm_phid_sq - ip_xphid^2 ;
            b = 2*( normx_minus_eps(l)*ip_phih_phid  -  ip_xphiH(l)*ip_xphid ) ;
       
            term1 = normx_minus_eps(l)*ip_phih_phid^2 ;
            term2 = 2*ip_xphid*ip_xphiH(l)*ip_phih_phid ;
            term3 = norm_phiH_sq(l)*ip_xphid^2 ;

            c = term1 - term2 + term3 ;
            c = c/norm_phid_sq ;

           
            % computing roots
            root1 = -b + sqrt(b^2 - 4*a*c) ;
            root1 = root1/(2*a) ;

            ip_phid_phihgamma_1 = ip_phih_phid + root1*norm_phid_sq ;
            ip_xphihgamma_1     = ip_xphiH(l) + root1*ip_xphid ;

            root1_check = (normx_minus_eps(l)*ip_phid_phihgamma_1/ip_xphid) - ip_xphihgamma_1 ;
        
            root2 = -b - sqrt(b^2 - 4*a*c) ;
            root2 = root2/(2*a) ;

            % to select the correct root
            if root1*root2 >= 0
                if (root1_check >= 0)
                    step_size(l) = root1 ;
                else
                    step_size(l) = root2 ;
                end
            else
                step_size(l) = max(root1 , root2) ;
            end
   
        end

%% Update

        H(:,l)            = H(:,l) + step_size(l)*( G_oracle_out(:,l) - H(:,l) ) ;

        ip_xphiH(l)       = ip_xphiH(l) + step_size(l)*(ip_xphig - ip_xphiH(l)) ;

        phiadjphi_H(:,l)  = phiadjphi_H(:,l) + step_size(l)*( phiadjphi_G(:,l) - phiadjphi_H(:,l) ) ;

        norm_phiH_sq(l)   = norm_phiH_sq(l) + step_size(l)^2*norm_phid_sq + 2*step_size(l)*ip_phih_phid ;

    end

    %% collecting data of the Plot-signal for plotting purpose

    % stepsize
    gamma(t) = step_size(Plot_signal); 

    % eta value for sub-optimality
    eta(t) = eta_val(Plot_signal) ;

    % first-order sub-optimality
    optimality_check(t) = 0.5 + (eta_grad(:,Plot_signal)'*H(:,Plot_signal))/(2*tau*norm(eta_grad(:,Plot_signal),1)) ;

    % iterate h
    F(:,t) = eta_val(Plot_signal).*H(:,Plot_signal) ;


end

end


%% Different descent direction Oracles in FLIPS

% For this part, the function ProjectionOntoL1Ball can be downloaded from:
% John Duchi (jduchi@cs.berkeley.edu) https://stanford.edu/~jduchi/projects/DuchiShSiCh08.html


% Simple Quadratic Oracle
function [g] = Simple_quad_oracle(h, grad, betainv,tau)

% gradient step for any LIP
g = h - betainv*grad ;

% LIP specific oracle, for sparse regression problems where c = l1-norm

g = l_inf_projection_oracle(g, tau); % computes projection onto l1-balls of radius tau

end


% Accelerated Quadratic Oracle
function [g, d_new] = Accelerated_quad_oracle(h, grad, d_old, betainv, momentum_para,tau) 

% gradient step for any LIP
g = h - betainv*grad - (betainv*momentum_para)*d_old ;

% LIP specific oracle, for sparse regression problems where c = l1-norm

g = l_inf_projection_oracle(g, tau) ; % computes projection onto l1-balls of radius tau
d_new = g - h ;

end


% Projection onto l1-ball

function w = l_inf_projection_oracle(v, tau)
% w = ProjectOntoL1Ball(v, b) returns the vector w which is the solution
%   to the following constrained minimization problem:
%
%    min   ||w - v||_2
%    s.t.  ||w||_inf <= tau.
%
%   It performs Euclidean projection of v to the l_inf-ball of radius tau.

if (tau < 0)
  error('Radius of L_inf ball is negative: %2.3f\n', tau);
elseif (norm(v, inf) <= tau)
  w = v;
else
  w = sign(v).*min(abs(v), tau) ;
end

end


