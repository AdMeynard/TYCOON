function crit = functional_value(F,sig,Fs,mu,nu,alpha,gamma)
% Estimation of the functional value
% Entries:
%   - F: the Ideal TF estimation
%   - sig: time samples of the signal
%   - Fs: Sampling Frequency
%   - mu: mu hyperparameter
%   - nu: nu hyperparameter
%   - alpha: alpha parameter
%   - gamma: gamma hyperparameter
%
N = length(sig);
M = ceil(N/2);


omega = linspace(0,Fs/2,M);% Frequency axis
xi = linspace(-Fs/2,Fs/2,N);% frequency axis for the FFT of F wrt time

crit1 = (1/Fs)*sum( (real((Fs/2/(M-1))*sum(F)) - sig(:).').^2);

if norm(alpha) == 0
    fftT_F = fftshift(fft(F.'),1).'; %  FFT wrt t for F(t,omega)
    crit2 = 1/(2*(M-1)*N) * sum(sum(4*pi^2*abs((repmat(xi,M,1) - repmat(omega',1,N)).*fftT_F).^2));
else
    crit2 = 1/(2*(M-1)) * sum(sum( abs( deriv_t(F,Fs,N,M) - 2i*pi*bsxfun(@times,omega',F) + bsxfun(@times,alpha,deriv_f(F,Fs,N,M)) ).^2 ));
end

crit3 = regul_Lp_value(F,1);

crit = crit1 + mu*crit2 + nu*crit3 + 0.5*gamma*norm(alpha(:)).^2;
end
