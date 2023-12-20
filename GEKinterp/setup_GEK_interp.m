function [v_GEK, beta] = setup_GEK_interp(X, Y, theta)
%
% compute the fundamental quantities for setting up a GEK predictor
%
% Inputs
% X     = sample points
% Y     = sample values
% theta = hyper parameters: correlation length
%
% Outputs
% v_GEK = GEK vector
% beta  = regression coeff (constant regression only)
%-----------------------------------------------------

dim   = size(X,1);
n     = size(X,2);
% make sure to use sample locations as column vectors


% set user parameters
use_grads  = 1;
regularize = 0;
corrmodel  = 2; %  1 for GAUSS, 2 for CUBIC

% set up correlation matrix
[R] = corr_matrix(X, theta, n, dim, use_grads, regularize, corrmodel);

% compute predictor data
F  = [ones(n,1);zeros(n*dim,1)];
v1 = linsolve(R,Y);
v2 = linsolve(R,F);

beta = (F'*v1)/(F'*v2);

v_GEK = v1 - beta*v2;

return;
end