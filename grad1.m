function grad1 = grad1(F,s,Fe)
% Compute the first part of the gradient

    N = length(s);
    M = ceil(N/2);
    grad1 = (2*pi/M)*repmat((pi*Fe/M)*sum(real(F)) - s(:).',M,1);
    
end

