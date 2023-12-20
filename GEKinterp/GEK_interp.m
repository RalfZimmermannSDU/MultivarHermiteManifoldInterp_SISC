function [yinterp] = GEK_interp(xstar, X, v_GEK, beta, theta, use_grads, corrmodel)
%------------------------------------------------------------------------------
%
% evaluate the gradient-enhanced Kriging predictor
%  
% INPUTS:
%   xstar     : trial point 
%   X         : matrix of sample points
%   v_GEK     : Kriging vector
%   beta      : regression coeff
%   theta     : distants weights
%   use_grads : boolean: gradients available? 0/1
%   corrmodel : Which correlation model?
%
% OUTPUTS:
%   yinterp   : interpolated output at xstar
%------------------------------------------------------------------------------

dim   = size(X,1);
n     = size(X,2);

if corrmodel ==1
            rxstar = GAUSS_corr_vector(X, xstar, theta, n, dim, use_grads);
else
            rxstar = CUBIC_corr_vector(X, xstar, theta, n, dim, use_grads);
end
yinterp = beta + v_GEK'*rxstar;

end
