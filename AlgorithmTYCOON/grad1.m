function grad1 = grad1(F,s,Fs)
%GRAD1 Compute the first part of the gradient
% Entries:
%   - F: the Ideal TF estimation
%   - s: time samples of the signal
%   - Fs: Sampling Frequency

    N = length(s);
    M = ceil(N/2);
    grad1 = 1/(M-1) * repmat( (Fs/2/(M-1))*sum(real(F)) - s(:).' , M, 1);
    
end

