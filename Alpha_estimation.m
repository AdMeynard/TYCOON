function alpha = Alpha_estimation(F,Fs,h,mu,G)
% Alpha estimation
% Entries:
%   - F: the TF matrix
%   - Fs: sampling frequency
%   - h: regularization value
%   - mu: mu hyperparameter
%   - G: 1 if G estimation, alpha estimation otherwise

[M,N] = size(F);


omega = linspace(0,Fs/2,M);% Frequency axis
D1 = deriv_t(F, Fs, N, M ) - 1i*2*pi*bsxfun(@times,omega',F);
D2 = deriv_f(F,Fs,N,M);

if G==1
    alpha = -mu*(D2'*D1)./(mu*D2+h);
  else   
    d =  diag(real(D2'*D1));
    d = -mu*d(:)';
    alpha = d./(mu*sum(abs(D2).^2,1)+h);
end

end