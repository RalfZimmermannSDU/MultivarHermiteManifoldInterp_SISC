function [ ] = FD_test_coeffs(Wlocs, ...
                              d1coeffs,...
                              d2coeffs,...
                              delta,...
                              theta,...
                              N,...
                              dim,...
                              index)
% Finite difference (FD) check of Kriging interpolator:
% Are the derivatives of the coefficient functions reproduced?
% INPUT
% Wlocs : sample locations in parameter space
% delta : FD step size
% theta : weights for GKE (= gradient enhanced Kriging)
% N     : number of sample points
% dim   : dimension of paramter space
% index : check and interpolate the sample point "index"

disp(['Perform FD check for interpolation coefficient functions'])
disp(['at sample location ', num2str(index), ':'])

% set up interpolation function 
Yphi            = zeros(N*(dim+1),1);
Yphi(index)     = 1.0;
Yphi(N+1:2*N)   = d1coeffs(index,:);
Yphi(2*N+1:3*N) = d2coeffs(index,:);
% set up predictor data
[v_GEK, beta] = setup_GEK_interp(Wlocs, Yphi, theta);

wstar = Wlocs(:,index);  % this is the sample location
                         % at which the partial derivatives are checked
% evaluate interpolator at wstar
[yinterp] = GEK_interp(wstar, Wlocs, v_GEK, beta, theta, 1, 2)


% test if the derivatives are reproduced
% A) coordinate direction 1
d1ystar = zeros(N,1);
for i=1:N
    xstar   = Wlocs(:,i);
    xstar_p = xstar + delta*[1;0];   % for central FD
    xstar_m = xstar - delta*[1;0];

    rxstar_p = CUBIC_corr_vector(Wlocs, xstar_p, theta, N, 2, 1);
    rxstar_m = CUBIC_corr_vector(Wlocs, xstar_m, theta, N, 2, 1);

    yinter_p = beta + v_GEK'*rxstar_p;
    yinter_m = beta + v_GEK'*rxstar_m;
    d1ystar(i) = (1.0/(2*delta)) * (yinter_p - yinter_m);
end
% B) coordinate direction 2
d2ystar = zeros(N,1);
for i=1:N
    xstar   = Wlocs(:,i);
    xstar_p = xstar + delta*[0;1];   % for central FD
    xstar_m = xstar - delta*[0;1];

    rxstar_p = CUBIC_corr_vector(Wlocs, xstar_p, theta, N, 2, 1);
    rxstar_m = CUBIC_corr_vector(Wlocs, xstar_m, theta, N, 2, 1);

    yinter_p = beta + v_GEK'*rxstar_p;
    yinter_m = beta + v_GEK'*rxstar_m;
    d2ystar(i) = (1.0/(2*delta)) * (yinter_p - yinter_m);
end

sampled_vs_interpolated = [Yphi(N+1:2*N), d1ystar, Yphi(2*N+1:3*N), d2ystar]

% END: test if the derivatives are reproduced

disp('END: Perform FD check for coefficient functions.')

%
end