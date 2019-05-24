function crit = regul_Lp_value(F,ord)
%REGUL_LP_VALUE Compute ||F||_p^p = sum_{t,f}^{N,M} |F_{t,f}|^p
% usage:	crit = regul_Lp_value(F,ord)

% Entries:
%   - F: the TF matrix
%   - ord: value of p (1 or 2)
%

crit=norm(F(:),ord)^ord;

if 0
    [nFreq, nTime] = size(F) ;
    tmp = zeros(1, nTime) ;
    
    for ii = 1: nTime
        tmp(ii) = norm(F(:, ii), ord).^ord ;
    end
    
    crit = norm(tmp,2) ; % L^2 norm in the time axis
end
end

