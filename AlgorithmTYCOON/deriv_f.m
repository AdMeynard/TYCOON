function [ Df_F ] = deriv_f(F, Fs, N, M )
%DERIV_F derivation with respect to frequency using Fourier transform
% Entries:
%   - F: the TF matrix
%   - Fs: sampling frequency
%   - N: length of the time vector
%   - M: length of the frequency vector


Df_F = deriv_t(F.', 2*(M-1)/Fs,M, N).';

end

