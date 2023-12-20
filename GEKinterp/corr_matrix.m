%------------------------------------------------------------------------------
% control function for setting up correlation matrix
% only rudimentary "cubic" and "Gauss" correlation kernels are hard-coded
%------------------------------------------------------------------------------
function [R] = corr_matrix(X, theta, n, dim, use_grads, regularize, corrmodel)

    if corrmodel == 1
        % use 'GAUSS':
        R = GAUSS_corr_matrix(X, theta, n, dim, use_grads, regularize);
    elseif corrmodel == 2
        % use 'CUBIC':
        R = CUBIC_corr_matrix(X, theta, n, dim, use_grads, regularize);
    else
        disp('WARNING: correlation model not available!')
        disp('=> Use Gaussian instead.')
        corrmodel = 1;
        R = GAUSS_corr_matrix(X, theta, n, dim, use_grads, regularize);
    end
return;
%------------------------------------------------------------------------------
end