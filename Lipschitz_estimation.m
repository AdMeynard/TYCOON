function [lambda,x] = Lipschitz_estimation(s,Fs,mu,alpha,Lip0,x0,stopping_threshold)
% Estimation of the Lipschitz constant using power iteration
% Entries:
%  - s: time samples of the signal
%  - Fs: Sampling frequency
%  - mu: the mu hyper-parameter
%  - alpha: value of alpha parameter
%  - Lip0: initialization for the Lipschitz constant
%  - x0: vector initialization
%  - stopping_threshold: to check convergence of the algorithm

N = length(s);
M = ceil(N/2);
omega = linspace(0,Fs/2,M); % Frequency axis

x = x0;
lambdaprec = Lip0+1;
lambda = Lip0;

while (abs(lambda-lambdaprec)/lambda>stopping_threshold) %stopping criterion
    lambdaprec = lambda;
    z = x/lambda; %normalization
    x = grad1(z,s,Fs) + mu*grad2(z,N,omega,M,Fs,alpha); % Linear Operator which the biggest eigen value is sought
    lambda = norm(x(:),'inf'); % Inf norm of the biggest eigen value
    %disp(abs(lambda-lambdaprec)/lambda)
end

end

