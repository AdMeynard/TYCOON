function [ Dt_F ] = deriv_t(F, Fs, N, M )

xi = linspace(-Fs/2,Fs/2,N);
fftT_F = fftshift(fft(F.'),1).'; % FFT wrt t of F(t,omega)


fft_Dt_F = 2*pi*1i*(repmat(xi,M,1)).*fftT_F;

Dt_F = ifft(ifftshift(fft_Dt_F.',1)).';


end

