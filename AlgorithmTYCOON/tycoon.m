function [F,alpha] = tycoon(s,Fs,mu,nu,alpha_regul,F0,alpha0,stop_eps, DEBUG)
%TYCOON Time-frequency bY Convex OptimizatiON algorithm
%   Usage: 
%    [F,alpha] = tycoon(s,Fs,mu,nu,alpha_regul,F0,alpha0,stop_eps, DEBUG)
%   Entries:
%   - s: time samples of the signal
%   - Fs: Sampling Frequency
%   - mu: the mu hyperparameter
%   - nu: the nu hyperparameter
%   - alpha_regul: the hyperparameter on ||alpha||
%   - F_0: Initial guess of the Ideal TF estimation
%   - alpha0: Initial guess of alpha
%   - stop_eps: stopping criterion
%   - DEBUG: if DEBUG=1, display of the stopping criterion through iterarions
% 
%   Output:
%   - F: estimated TF matrix
%   - alpha: estimated chirp factor


[M,N] = size(F0);
Lipschitz_constant = 1;
eiv = ones(M,N);
F=F0;
alpha=alpha0;
G = 0;
ok = 1;
k = 0;

while ok
    k=k+1;
    %descent on  F
    Fold = F;
    alphaold = alpha;
    [Lipschitz_constant,eiv] = Lipschitz_estimation(s,Fs,mu,alpha,Lipschitz_constant,eiv,1e-4);% Restimate the Lipschitz constant
    [F] = proximal_descent(s,Fs,F,mu,1.1*Lipschitz_constant,nu,alpha,alpha_regul,stop_eps, DEBUG);
    
    
    %descent on alpha
    alpha = Alpha_estimation(F,Fs,alpha_regul,mu,G);
    
    
    stop_critF = norm(F(:)-Fold(:))/norm(F(:));
    stop_critAlpha = norm(alpha(:)-alphaold(:))/norm(alpha);
    if DEBUG
        fprintf('[In TYCOON] it: %d, stop_critF = %f, stop_critAlpha = %f\n',k,stop_critF,stop_critAlpha);
    end
    
    if ((norm(F(:)) < 1e-4) || ((stop_critF < stop_eps) && stop_critAlpha < stop_eps))
        ok = 0;
    end
    
end

