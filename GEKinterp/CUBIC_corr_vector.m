
function [rx] =  CUBIC_corr_vector(X, xstar, theta, n, dim, use_grads)
%------------------------------------------------------------------------------
%
%CUBIC_corr_vector:
%    Set up correlation vector w.r.t CUBIC correlation kernel
% INPUTS:
%   X         = matrix of sample points
%   xstar     = trial point. 
%   theta     = distants weights
%   n         = number of sample points
%   dim       = dimension of sample point
%   use_grads = boolean: gradients available? 0/1
%
% OUTPUTS
%  rx         = correlation vector r(xstar)
%------------------------------------------------------------------------------
    % initialize corr vector rx
    rx = zeros(n*(use_grads*dim+1),1);

    % compute inter-correlations block
    for i=1:n % in range(n):
        rx(i) = CUBIC_corr(theta, X(:,i), xstar, dim);
    end

    if use_grads
        % compute cross-correlations block
        for k=1:dim % in range(dim):
            for i=1:n % in range(n):
                rx(i+k*n) = d_k_CUBIC_corr(k, theta, X(:,i), xstar, dim);
            end
            %end i
        end
        %end k
    end

    return;
end
% END: CUBIC_corr_vector
%------------------------------------------------------------------------------