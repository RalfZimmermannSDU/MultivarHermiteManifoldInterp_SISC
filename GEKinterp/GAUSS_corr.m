function [r] = GAUSS_corr(theta, xi, xj, dim)
%------------------------------------------------------------------------------
% compute GAUSS correlation between vectors xi, xj of dimension dim
%------------------------------------------------------------------------------

r = 0.0;
for k = 1:dim
    r = r + theta(k) * (xi(k)-xj(k))*(xi(k)-xj(k));
end

r = exp(-r);
return;
%------------------------------------------------------------------------------
end