
function [R] =  GAUSS_corr_matrix(X, theta, n, dim, use_grads, regularize)
%------------------------------------------------------------------------------
%GAUSS_corr_matrix Set up correlation matrix w.r.t GAUSS
% INPUTS:
%   X         : matrix of sample points
%   theta     : distants weights
%   n         : number of sample points
%   dim       : dimension of sample point
%   use_grads : boolean: gradients available? 0/1
%   regularize: boolean:0/1, add small epsilon to diagonal?
%
% OUTPUTS:
%    R        : Gauss correlation matrix
%------------------------------------------------------------------------------
    % initialize corr matrix R
    if use_grads
        R = zeros(n*(dim+1), n*(dim+1));
    else
        R = zeros(n,n);
    end

    % compute auto-correlations block
    for i =1:n
        for j=1:n
            R(i,j) = GAUSS_corr(theta, X(:,i), X(:,j), dim);
        end
    end

    if use_grads
        % i.e.   |   R      d(:,1)R   d(:,2)R  |
        %     R =|d(1,:)R  d^2(1,1)R d^2(1,2)R |
        %        |d(2,:)R  d^2(2,1)R d^2(2,2)R |
        
        %compute cross-correlations blocks
        % d(k,:)R(xi, xj), k = 1,...,dim
        for k=1:dim % in range(dim):
            for i=1:n % in range(n):
                for j=1:n % in range(n):
                    R(i+k*n,j) = R(i,j)*(-2)*theta(k)*(X(k,i)-X(k,j));
                end
            end
        end

        
        %compute cross-correlations blocks
        % d(:,k)R(xi, xj), k = 1,...,dim
        for k=1:dim % in range(dim):
            for i=1:n % in range(n):
                for j=1:n % in range(n):
                    R(i,j+k*n) = R(i,j)*(-2)*theta(k)*(X(k,j)-X(k,i));
                end
            end
        end
  

        % compute auto-correlations between derivatives
        % d(l,k)R(xi,xj), k,l = 1,..., dim
        for l=1:dim % in range(dim):
            for k=1:dim % in range(dim):
                for i=1:n % in range(n):
                    for j=1:n % in range(n):
                        if (l~=k)
                            R(i+l*n , j+k*n) = ...
                                R(i,j)*(-4)*theta(l)*theta(k)*(X(l,i)-X(l,j))*(X(k,i)-X(k,j));
                        else
                            R(i+l*n, j+k*n) = ...
                                R(i,j)*2*theta(k)*(1-2*theta(k)*(X(k,i)-X(k,j))*(X(k,i)-X(k,j)));
                        end
                        %End "if l~=k"
                    end
                    %End "for j..."
                end
                %End "for i..."
            end
            %End "for k..."
        end
        %End "for l..."
    end
    %End "if use_grads"

    % Regularization
    if regularize
        epsilon = 1.0e-8;
        R = R + epsilon*eye(size(R));
    end
return;
end
%------------------------------------------------------------------------------