function [F] = proximal_descent(sig,Fs,F_0,mu,Lipschitz_constant,nu,alpha,gamma,stop_eps,DEBUG)
%PROXIMAL_DESCENT Functional optimization using proximal descent algorithm
% Entries:
%   - sig: time samples of the signal
%   - Fs: Sampling Frequency
%   - F_0: Initial guess of the Ideal TF estimation
%   - mu: the mu hyperparameter
%   - Lipschitz_constant: value of the Lipschitz constant
%   - nu: the nu hyperparameter
%   - alpha: value of alpha
%   - gamma: the gamma hyperparameter
%   - stop_eps: stopping criterion
%   - DEBUG: if DEBUG=1, display of the stopping criterion through iterations


N = length(sig);
M = ceil(N/2);
omega = linspace(0,Fs/2,M); % Frequency axis

F = F_0; % the TF matrix
z = F_0;

ok = 1;
k=0;

while ok
    k = k+1;
    
    crit_now = functional_value(F,sig,Fs,mu,nu,alpha,gamma);
    
    F_old = F;
    gradz = grad1(z,sig,Fs)+mu*grad2(z,N,omega,M,Fs,alpha); % Gradient

    for tau=linspace(Lipschitz_constant/100,Lipschitz_constant,10) % step estimation
        FF = z - (1/tau)*gradz;
        FF = prox(FF,nu,tau,1); % Next iteration
        newcrit = functional_value(FF,sig,Fs,mu,nu,alpha,gamma);
        if newcrit <= crit_now
            break;
        end
    end
    if newcrit > crit_now % Monotonic version
        F = F_old;
    else
        F = FF;
    end
    
    gamma1 = (k+1)/(k+2);
    gamma2 = (k-1)/(k+3);
    z = F + gamma1*(FF-F) + gamma2*(F - F_old); % relaxation (FISTA)
    
    crit_stop = norm(F(:)-F_old(:))/norm(F(:));
    if DEBUG
       fprintf('  (in prox descent) it: %d, functional_value = %f, stopping criterion: %f\n',k,crit_now,crit_stop);
    end
    
    if FF == 0
        break ;
    end
    if crit_stop < stop_eps
        ok = 0;
    end
    
end



