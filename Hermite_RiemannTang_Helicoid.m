%
%
% This script performs the experiment from Section 5.1 of the paper
% "MULTIVARIATE HERMITE INTERPOLATION ON RIEMANNIAN MANIFOLDS",
%  Zimmermann/Bergmann
%
%
clear;
%
%
%
%
addpath('GEKinterp/')
addpath('aux4test/')
%
%
%*** USER PARAMETERS ***
n1    = 3   % sample resolution in coordinate direction 1
n2    = 3   % sample resolution in coordinate direction 1
res   = 101 % resolution for trial space
delta = 1.0e-5     % delta for finite difference tests
tau   = 1.0e-8     % convergence threshold for Riemannian optimization
index = 3          % check interpolation condition for sample point "index"
theta = [0.5,0.5]; % correlation weights for Kriging
%
%
% create sample points 
%
%
%
a = -0.25*pi;
b =  0.25*pi;
%
%
wspace1  = linspace(a,b,n1); 
wspace2  = linspace(a,b,n2);
[W1, W2] = meshgrid(wspace1, wspace1);
Wspace   = [W1(:)'; W2(:)'];

% trial space
wspace1_trial = linspace(a,b,res);
wspace2_trial = linspace(a,b,res);

dim = 2; % 2D sample space

N = n1*n2;

plot_helicoid(wspace1, wspace2);

set(0,'defaultTextInterpreter','latex'); %trying to set the default
hold on

% choose test function
testfun_handle = @(x,y) testfun_gauss_S2(x,y);

% create all data required for Hermite interpolation
tic;
[Wlocs, samplocs, d1samplocs, d2samplocs, pstar, logs_array, d1_logs_array, d2_logs_array] = create_sample_data_Sphere_tang(testfun_handle, Wspace, N, dim);
% dimension of matrices
d = size(samplocs, 1)

% set up sample data vecs for interpolation
Yphi = eye(N*(dim+1),N*(dim+1));
% each column j of Yphi contains
% * the sample values of the coeff funcs
% * the partial derivatives in direction d1 of the coeff funcs
% * the partial derivatives in direction d2 of the coeff funcs
% for the linear combination
%
% yhat = \sum_j \lambda_j(w) logs_array(j) + \sum_ij \mu_ij(w) di_logs_array(j)
%
% the sample data vector is one column of Yphi, namely
% Yphi(:, k) = (lambda_k(w1),...,lambda_k(wn),
% d1lambda_k(w1),...,d1lambda_k(wn), ...)

% Interpolate the coefficient functions 
% and solve weighted barycenter problem

% set up Kriging predictor for each Yphi(:,k)
v_GEK_array = zeros(N*(dim+1), N*(dim+1));
beta_array  = zeros(N*(dim+1),1);
for k=1:N*(dim+1)
    [v_GEK, beta] = setup_GEK_interp(Wlocs, Yphi(:,k), theta);
    v_GEK_array(:,k) = v_GEK;
    beta_array(k)    = beta;
end
t_prep = toc;
disp(['data preprocessing: ', num2str(t_prep), 's'])


yinterp = zeros(d,d,res,res);

error_field = zeros(res,res);

%
disp("Running interpolation over trial space, please wait.")
tic;
for j=1:res
    for k=1:res
        wstar = [wspace1_trial(j); wspace2_trial(k)];
        
        % true function value:
        [q_true, d2q_true, d2q_true] = testfun_handle(wspace1_trial(j), wspace2_trial(k));
        
        % interpolate the coefficient funcs
        weights = zeros(N*(dim+1),1);
        for l=1:N*(dim+1)
            weights(l) = GEK_interp(wstar, Wlocs, v_GEK_array(:,l), beta_array(l), theta, 1, 2);
        end
        % Do the weight functions sum up to 1?
        %disp(['sum of sample weights = ', num2str(sum(weights(1:N)))])
        
        % interpolate
        y_tang = zeros(d,1);
        for l=1:N
            y_tang = y_tang + weights(l)*logs_array(:,l) + weights(N+l)*d1_logs_array(:,l) + ...
                     weights(2*N+l)*d2_logs_array(:,l);
        end
        % map y_tang to manifold
        yinterp(:,j,k) = exp_sphere(pstar,1,  y_tang);

        % record error
        error_field(j,k) = norm(yinterp(:,j,k) - q_true, 2);
    end
end
t_interp = toc;
disp(['data interpolation: ', num2str(t_interp), 's'])


disp([' ']);
disp(['*****************************************']);
disp([' Av. interpolation error             : ', num2str(sum(sum(error_field))/(res*res))])
disp([' Max interpolation error             : ', num2str(max(max(error_field)))])
disp(['*****************************************']);


% plot errors:
[Xerrors, Yerrors] = meshgrid(wspace1_trial, wspace2_trial);

subplot(1,2,2)
surf(Xerrors, Yerrors, error_field)
xlabel('$\omega_1$','interpreter','latex')
ylabel('$\omega_2$','interpreter','latex')
zlabel('Errors $\|f(\omega) - \hat f(\omega)\|_2$','interpreter','latex')
colormap gray;
title('Interpolation errors: TSHI', 'interpreter','latex')
brighten(0.8)

%
%
%
%



%***************************************
disp([' ']);
disp(['**************************************************']);
disp('Perform FD check for interpolated manifold function')
disp(['**************************************************']);

d1_error = zeros(n1*n2,1);
d2_error = zeros(n1*n2,1);
for j=1:n1*n2
    index = j;

    wstar = Wlocs(:,index);
    wstar_1p = wstar + delta*[1;0];
    wstar_1m = wstar - delta*[1;0];
    wstar_2p = wstar + delta*[0;1];
    wstar_2m = wstar - delta*[0;1];
    
    % interpolate the coefficient funcs
    weights = zeros(N*(dim+1),1);
    for l=1:N*(dim+1)
        weights_1p(l) = GEK_interp(wstar_1p, Wlocs, v_GEK_array(:,l), beta_array(l), theta, 1, 2);
        weights_1m(l) = GEK_interp(wstar_1m, Wlocs, v_GEK_array(:,l), beta_array(l), theta, 1, 2);
        weights_2p(l) = GEK_interp(wstar_2p, Wlocs, v_GEK_array(:,l), beta_array(l), theta, 1, 2);
        weights_2m(l) = GEK_interp(wstar_2m, Wlocs, v_GEK_array(:,l), beta_array(l), theta, 1, 2);
    end
    %
    y_tang_1p = zeros(d,1);
    y_tang_1m = zeros(d,1);
    y_tang_2p = zeros(d,1);
    y_tang_2m = zeros(d,1);
    for l=1:N
        y_tang_1p = y_tang_1p + weights_1p(l)*logs_array(:,l) + weights_1p(N+l)*d1_logs_array(:,l) + ...
                    weights_1p(2*N+l)*d2_logs_array(:,l);
        y_tang_1m = y_tang_1m + weights_1m(l)*logs_array(:,l) + weights_1m(N+l)*d1_logs_array(:,l) + ...
                    weights_1m(2*N+l)*d2_logs_array(:,l); 
        y_tang_2p = y_tang_2p + weights_2p(l)*logs_array(:,l) + weights_2p(N+l)*d1_logs_array(:,l) + ...
                    weights_2p(2*N+l)*d2_logs_array(:,l);
        y_tang_2m = y_tang_2m + weights_2m(l)*logs_array(:,l) + weights_2m(N+l)*d1_logs_array(:,l) + ...
                    weights_2m(2*N+l)*d2_logs_array(:,l);
    end

    % map to manifold
    y_1p = exp_sphere(pstar, 1, y_tang_1p);
    y_1m = exp_sphere(pstar, 1, y_tang_1m);
    y_2p = exp_sphere(pstar, 1, y_tang_2p);
    y_2m = exp_sphere(pstar, 1, y_tang_2m);

    d1ystar = (1.0/(2*delta)) *(y_1p - y_1m);
    d2ystar = (1.0/(2*delta)) *(y_2p - y_2m);

    d1_error(j) = norm(d1ystar - d1samplocs(:,index), 'fro')/norm(d1samplocs(:,index), 'fro');
    d2_error(j) = norm(d2ystar - d2samplocs(:,index), 'fro')/norm(d1samplocs(:,index), 'fro');
    %disp(['Error( d/d1 f(wj) - d/d1 hat f(wj)) at index j =', num2str(index), ': ',...
    %    num2str(d1_error(j))])

    %disp(['Error( d/d2 f(wj) - d/d2 hat f(wj)) at index j =', num2str(index), ': ',...
    %    num2str(d2_error(j))])
end
av_FD_err1 = sum(d1_error)/(n1*n2)
av_FD_err2 = sum(d2_error)/(n1*n2)
