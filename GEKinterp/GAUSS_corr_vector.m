function [rx] = GAUSS_corr_vector(X, xstar, theta, n, dim, use_grads)
%------------------------------------------------------------------------------
%
%GAUSS_corr_vector:
%    Set up correlation vector w.r.t GAUSS correlation kernel
% INPUTS:
%   X         : matrix of sample points
%   xstar     : trial point. 
%   theta     : distants weights
%   n         : number of sample points
%   dim       : dimension of sample point
%   use_grads : boolean: gradients available? 0/1
%   regularize: boolean:0/1, add small epsilon to diagonal?
%
% OUTPUTS:
%    rx       : Gauss correlation vector r(x)
%------------------------------------------------------------------------------
    % initialize corr vector rx
    rx = zeros(n*(use_grads*dim+1),1);

    % compute inter-correlations block
    for i=1:n % in range(n):
        rx(i) = GAUSS_corr(theta, X(:,i), xstar, dim);
    end
    if use_grads
        % compute cross-correlations block
        for k=1:dim % in range(dim):
            for i=1:n % in range(n):
                rx(i+k*n) = rx(i)*(-2)*theta(k)*(X(k,i)-xstar(k));
            end
            %end i
        end
        %end k
    end

return;
end
% END: GAUSS_corr_vector
%------------------------------------------------------------------------------