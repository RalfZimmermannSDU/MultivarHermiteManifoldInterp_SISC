% *****************************************************************************
% Numerical optimzation for computing Riemann centers on sphere
% simple gradient descent
% Inputs: 
%    Locs = list locations on manifold
%           2D array of dims [3,n] n= number of samples, 3=embedding
%           dimension
% weights = vector of weights, same length das Locs list
%     tau = convergence threshold
%
% Output:
%   Vstar = Riemannian center, solving 
%           min (sum weights[j] dist(Uj, Vstar)^2)
%
% *****************************************************************************
function [Vstar, count, fail] =  findCenter_Sphere(Locs, weights, q0, tau)
    %0. get dimensions
    [d,n] = size(Locs);
    
    %**************************************************************************
    %1. preprocessing: compute initial guess
    %**************************************************************************
    %compute the arithmetic mean in a tangent space of one of the data points
    nr_Locs = length(weights);

    %initialize gradient (actually, this is the negative of the gradient)
    grad = zeros(d,1);
    Vstar = zeros(d,1);
    for l=1:n
        Delta_l = log_sphere(q0, Locs(:,l));
        grad    = grad + weights(l)*Delta_l;
    end
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
        Vstar = exp_sphere(Vstar, 1.0, alpha*grad);
        count = count +1;
        %update gradient
        grad = zeros(d,1);
        for l=1:n 
            Delta_l = log_sphere(Vstar, Locs(:,l));
            grad    = grad + weights(l)*Delta_l;
        end
        norm_grad = norm(grad, 2);
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