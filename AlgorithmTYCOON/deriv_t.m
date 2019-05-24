function [ Dt_F ] = deriv_t(F, Fs, N, M )
%DERIV_F derivation with respect to time using Fourier transform
% Entries:
%   - F: the TF matrix
%   - Fs: sampling frequency
%   - N: length of the time vector
%   - M: length of the frequency vector

xi = (-N/2:N/2-1)*(Fs/N);
fftT_F = fftshift(fft(F.'),1).'; % FFT wrt t of F(t,omega)


fft_Dt_F = 2i*pi*(repmat(xi,M,1)).*fftT_F;

Dt_F = ifft(ifftshift(fft_Dt_F.',1)).';


end

