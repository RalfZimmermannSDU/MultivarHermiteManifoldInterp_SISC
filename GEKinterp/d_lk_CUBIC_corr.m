function [r] =  d_lk_CUBIC_corr(l, k, theta, xi, xj, dim)
%------------------------------------------------------------------------------
% compute second order derivative of cubic correlation d(l,k)r(theta, xi,xj)
% 
% INPUTS
%  l, k : compute partial derivative with respect to variables l, k
% theta : correlation weights
%  xi   : sample location x_i
%  xj   : sample location x_j
%  dim  : dimension
% OUTPUTS
%  partial derivative d(l,k)r(theta, xi,xj)
%------------------------------------------------------------------------------
    r = 1.0;

    if l==k
        for h=1:dim % in range(dim):
            if (theta(h)*abs(xi(h)-xj(h)) < 1.0)
                if h==l
                    r = r *6*theta(h)*theta(h)*(1-2*theta(h)*abs(xi(h)-xj(h)));
                else
                    r = r *(1 - 3*theta(h)*theta(h)*(xi(h)-xj(h))*(xi(h)-xj(h)) ...
                     + 2*theta(h)*theta(h)*theta(h) * abs(xi(h)-xj(h))*abs(xi(h)-xj(h)) *abs(xi(h)-xj(h)));
                end
            else
                r = 0.0;
            end
        end
    else
        for h=1:dim % in range(dim):
            if (theta(h)*abs(xi(h)-xj(h)) < 1.0)
                if h == l
                    r = r * 6*theta(h)*theta(h)*sign(xi(h)-xj(h))*(theta(h)*(xi(h)-xj(h))*(xi(h)-xj(h)) - abs(xi(h)-xj(h)));
                elseif h == k
                    r = (-1)*r * 6*theta(h)*theta(h)*sign(xi(h)-xj(h))*(theta(h)*(xi(h)-xj(h))*(xi(h)-xj(h)) - abs(xi(h)-xj(h)));
                else
                    r = r*(1 - 3*theta(h)*theta(h)*(xi(h)-xj(h))*(xi(h)-xj(h)) ...
                     + 2*theta(h)*theta(h)*theta(h) * abs(xi(h)-xj(h)) *abs(xi(h)-xj(h))*abs(xi(h)-xj(h)));
                end
            else
                r = 0.0;
            end
        end
    end

return;
end
%------------------------------------------------------------------------------