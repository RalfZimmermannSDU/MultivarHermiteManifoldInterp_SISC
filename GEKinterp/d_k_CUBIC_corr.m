
function [r] =  d_k_CUBIC_corr(k, theta, xi, xj, dim)
%------------------------------------------------------------------------------
% compute partial derivative of 
% CUBIC correlation between vextors xi, xj of dimension dim
% derivative of cubic correlation d(k,:)r(theta, xi,xj)
% 
% INPUTS
%     k : compute partial derivative with respect to component k
% theta : correlation weights
%  xi   : sample location x_i
%  xj   : sample location x_j
%  dim  : dimension
% OUTPUTS
%  partial derivative d(k,:)r(theta, xi, xj)
%------------------------------------------------------------------------------
    r = 1.0;
    for l=1:dim % in range(dim):
        if (theta(l)*abs(xi(l)-xj(l)) < 1.0)
            if l == k
                r = r * 6*theta(l)*theta(l)*sign(xi(l)-xj(l))*(theta(l)*(xi(l)-xj(l))*(xi(l)-xj(l)) - abs(xi(l)-xj(l)));
            else
                r = r*(1 - 3*theta(l)*theta(l)*(xi(l)-xj(l))*(xi(l)-xj(l)) ...
                    + 2*theta(l)*theta(l)*theta(l) * abs(xi(l)-xj(l))*abs(xi(l)-xj(l)) *abs(xi(l)-xj(l)) );
            end
        else
            r = 0.0;
        end
    end
return;
end
%------------------------------------------------------------------------------