function Dot = performIND(ifi,ami,F,omega)

[M,N] = size(F);
nbComp = size(ifi,1);

alpha = omega(2)- omega(1) ;
itvPS = zeros(M,N) ;

for kk = 1:N
    for ll = 1:nbComp
        if ~isnan(ifi(ll,kk)./alpha)
            itvPS(min(M,max(1,round(ifi(ll,kk)./alpha))), kk) = ami(ll,kk)^2 ; 
        end
    end
end

Dot = slicedOT(itvPS, abs(F).^2) * 100;