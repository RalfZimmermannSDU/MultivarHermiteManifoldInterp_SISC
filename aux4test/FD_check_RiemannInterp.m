function [d1ystar, d2ystar] = FD_check_RiemannInterp(samplocs, ...
                                                     Wlocs, ...
                                                     q0, ...
                                                     tau, ...
                                                     v_GEK_array, ...
                                                     beta_array, ...
                                                     theta, ...
                                                     N, ...
                                                     dim, ...
                                                     delta, ...
                                                     index, ...
                                                     example)
% Make FD check, whether or not Riemannian interpolator matches gradients.
% Are the manifold derivatives of the test function reproduced?
%
% !Hard-coded for two-dimensional parameter space!
%
% INPUTS
% samplocs     : Riemannian samples
% Wlocs        : sample locations in parameter space
% q0           : start point for Riemannian optimization
% tau          : threshold for Riemannian optimization
% v_GEK_array  : array of v_GEK vectors
% beta_array   : araya of beta regression constants
% theta        : weights for GKE
% N            : number of sample points
% dim          : dimension of paramter space
% delta        : FD step size
% index        : check and interpolate the sample point "index"
% example      : 2-string: 'SO' or 'S2'
%              : 'S2' is for unit sphere in R^3 
%                (Section 5.1 of associated SISC paper)
%              : 'SO' is for special orthogonal group 
%                (Section 5.2 of associated SISC paper)
%
% OUTPUTS
% d1ystar      : partial derivative by first argument
% d2ystar      : partial derivative by second argument
%

wstar    = Wlocs(:,index);        % current trial point
wstar_1p = wstar + delta*[1;0];   % one step ahead
wstar_1m = wstar - delta*[1;0];   % one step back for central FD
wstar_2p = wstar + delta*[0;1];   % one step ahead
wstar_2m = wstar - delta*[0;1];   % one step back for central FD

% interpolate the coefficient funcs
weights_1p = zeros(N,1);
weights_1m = zeros(N,1);
weights_2p = zeros(N,1);
weights_2m = zeros(N,1);
for l=1:N
    weights_1p(l) = GEK_interp(wstar_1p, Wlocs, v_GEK_array(:,l), beta_array(l), theta, 1, 2);
    weights_1m(l) = GEK_interp(wstar_1m, Wlocs, v_GEK_array(:,l), beta_array(l), theta, 1, 2);
    weights_2p(l) = GEK_interp(wstar_2p, Wlocs, v_GEK_array(:,l), beta_array(l), theta, 1, 2);
    weights_2m(l) = GEK_interp(wstar_2m, Wlocs, v_GEK_array(:,l), beta_array(l), theta, 1, 2);
end

% find weighted center by solving Riemann optimization problem
if example == 'S2'
    [ystar_1p, count, fail_1p] = findCenter_Sphere(samplocs, weights_1p, q0, tau);
    [ystar_1m, count, fail_1m] = findCenter_Sphere(samplocs, weights_1m, q0, tau);
    [ystar_2p, count, fail_2p] = findCenter_Sphere(samplocs, weights_2p, q0, tau);
    [ystar_2m, count, fail_2m] = findCenter_Sphere(samplocs, weights_2m, q0, tau);
elseif example == 'SO'
    [ystar_1p, count, fail_1p] = findCenter_SOn(samplocs, weights_1p, q0, tau);
    [ystar_1m, count, fail_1m] = findCenter_SOn(samplocs, weights_1m, q0, tau);
    [ystar_2p, count, fail_2p] = findCenter_SOn(samplocs, weights_2p, q0, tau);
    [ystar_2m, count, fail_2m] = findCenter_SOn(samplocs, weights_2m, q0, tau);    
else
    disp("In FD_check_RiemannInterp.m: Example not set correctly.")
    disp("Choose string <SO> or string <S2>.")
end

if (fail_1p + fail_1m + fail_2p + fail_2m == 0)
    disp(['All Riemann opt. processes converged in -FD_check_RiemannInterp-!'])
else
    disp('Failure in -findCenter- in -FD_check_RiemannInterp-')
end

% central difference approximation
d1ystar = (1.0/(2*delta)) * (ystar_1p - ystar_1m);
d2ystar = (1.0/(2*delta)) * (ystar_2p - ystar_2m);
end  % End function