function [ Df_F ] = deriv_f(F, Fs, N, M )

tau = linspace(0,N/Fs,M);
tau = tau(:);
fftf_F = fft(F); % FFT wrt omega of F(t,omega)


fft_Df_F = 2*pi*1i*(repmat(tau,1,N)).*fftf_F;

Df_F = ifft(fft_Df_F);


end

