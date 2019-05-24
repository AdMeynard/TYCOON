function F_out = prox(F_in,nu,L,p_Lp)
%PROX Compute the proximity operator of L2 or L1 norm
%   F_out = argmin_F 0.5*||F_in - F||_2^2 + nu/L ||F||_p^p
% with p = {1,2}
%
% Entries:
%   - F_in: the TF matrix estimation to apply the proximity operator
%   - nu: the hyperparameter
%   - L: the Lipschitz constant
%   - p_Lp: value of p (1 or 2)
%

    switch p_Lp
        case 1
            F_out = sign(F_in).*max(0,abs(F_in)-nu/L);
        case 2
            F_out = F_in/(1+(2*nu/L));
    end
end

