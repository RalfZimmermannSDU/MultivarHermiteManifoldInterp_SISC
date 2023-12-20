% *****************************************************************************
% Numerical optimzation for computing Riemann centers on SOn
% simple gradient descent
% Inputs: 
%    Locs = list locations on manifold, here matrices are expected
%           3D array of dims [d,d,n] n= number of samples, 
%           dxd=embedding dimension
% weights = vector of weights, same length as Locs list
%     tau = convergence threshold
%
% Output:
%   Vstar = Riemannian center, solving 
%           min (sum weights[j] dist(Uj, Vstar)^2)
%
% *****************************************************************************
function [Vstar, count, fail] =  findCenter_SOn(Locs, weights, q0, tau)
    %0. get dimensions
    [d,d2,n] = size(Locs);
    
    %**************************************************************************
    %1. preprocessing: compute initial guess
    %**************************************************************************
    %compute the arithmetic mean in a tangent space of one of the data points
    nr_Locs = length(weights);

    %initialize gradient (actually, this is the negative of the gradient)
    grad = zeros(d,d);
    Vstar = zeros(d,d);
    for l=1:n
        Delta_l = log_SOn(q0, Locs(:,:,l));
        grad = grad + weights(l)*Delta_l;
    end
    grad = 0.5*(grad - grad');
    norm_grad = norm(grad, 2);
    %disp(['-find center- initial grad norm: ', num2str(norm_grad), '. Target:', num2str(tau)])
    

    %**************************************************************************
    %2. Numerical optimization
    %**************************************************************************

    count = 0;
    countmax = 2000;
    Vstar = q0;
    % step size
    alpha = 1.0;
    while (norm_grad>tau) && (count < countmax)
        % make a geodesic step in direction of gradient
        Vstar = exp_SOn(Vstar, alpha*grad);
        count = count + 1;
        %update gradient
        grad = zeros(d,d,1);
        for l=1:n 
            Delta_l = log_SOn(Vstar, Locs(:,:,l));
            grad = grad + weights(l)*Delta_l;
        end
        grad = 0.5*(grad - grad');
        norm_grad = norm(grad, 'fro');
        %disp(["current grad norm:", num2str(norm_grad)]);
    end

    if count < countmax
        %disp(['Function -findCenter- converged in ', num2str(count), ' iterations'])
        fail = 0;
    else
        %disp(['Function -findCenter- did not converged'])
        fail = 1;
    end
return;
end
% *********************************************