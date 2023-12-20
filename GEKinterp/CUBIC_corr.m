%------------------------------------------------------------------------------
% compute CUBIC correlation between vextors xi, xj of dimension dim
%------------------------------------------------------------------------------
function [r] =  CUBIC_corr(theta, xi, xj, dim)
    r = 1.0;
    for k=1:dim % in range(dim):
        if theta(k)*abs(xi(k)-xj(k)) < 1.0
            r = r*(1 - 3*theta(k)*theta(k)*(xi(k)-xj(k))*(xi(k)-xj(k)) ...
                + 2*theta(k)*theta(k)*theta(k)*abs(xi(k)-xj(k))*abs(xi(k)-xj(k)) *abs(xi(k)-xj(k))  );
        else
            r = 0.0;
        end
    end
return;
end
%------------------------------------------------------------------------------