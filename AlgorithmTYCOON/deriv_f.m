function [ Df_F ] = deriv_f(F, Fs, N, M )

Df_F = deriv_t(F.', 2*(M-1)/Fs,M, N).';

% tau = linspace(0,N/Fs,M);
% tau = tau(:);
% fftf_F = fft(F); % FFT wrt omega of F(t,omega)
% 
% 
% fft_Df_F = 2i*pi*(repmat(tau,1,N)).*fftf_F;
% 
% Df_F = ifft(fft_Df_F);


end

