%
% This script performs the experiment from Section 5.1 of the paper
% "MULTIVARIATE HERMITE INTERPOLATION ON RIEMANNIAN MANIFOLDS",
%  Zimmermann/Bergmann
%
clear;
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
index = 1          % check interpolation condition for sample point "index"
theta = [0.5,0.5]; % correlation weights
%
%
% create sample points 
%
%
%
a =-0.25*pi;
b = 0.25*pi;
%
wspace1 = linspace(a,b,n1);
wspace2 = linspace(a,b,n2);
% trial space
wspace1_trial = linspace(a,b,res);
wspace2_trial = linspace(a,b,res);

dim = 2; % 2D sample space

N = n1*n2;

% choose test function
testfun_handle = @(x,y) testfun_gauss_S2(x,y);

% create all data required for Hermite interpolation
tic;
[Wlocs, samplocs, d1samplocs, d2samplocs, d1coeffs, d2coeffs] = create_sample_data_Sphere(testfun_handle, wspace1,wspace2, N, dim);


% compute some q as a candidate for the barycenter
mean = sum(samplocs,2)/N;
q0 = mean/norm(mean);   % map to Sphere

% set up sample data vecs for interpolation
Yphi = zeros(N*(dim+1),N);
% each column j of Yphi contains
% * the sample values 
% * the partial derivatives in direction d1 
% * the partial derivatives in direction d2
% for the coefficient function j 
% in the objective function "sum_j phi_j(w) dist(q,pj)"
for k = 1:N
    Yphi(k,k) = 1.0;   % the phi's interpolate the unit vectors
    Yphi(N+1:2*N,k)   = d1coeffs(k,:);
    Yphi(2*N+1:3*N,k) = d2coeffs(k,:);
end

% Interpolate the coefficient functions 
% and solve weighted barycenter problem

% set up Kriging predictor for each Yphi(:,k)
v_GEK_array = zeros(N*(dim+1), N);
beta_array  = zeros(N,1);
for k=1:N
    [v_GEK, beta] = setup_GEK_interp(Wlocs, Yphi(:,k), theta);
    v_GEK_array(:,k) = v_GEK;
    beta_array(k)    = beta;
end
t_preprocessing = toc;
disp(['data preprocessing: ', num2str(t_preprocessing), 's'])

% plot data
% this opens a figure
plot_helicoid(wspace1, wspace2);

set(0,'defaultTextInterpreter','latex'); %trying to set the default
hold on

% FD check for coefficients
FD_test_coeffs(Wlocs, d1coeffs, d2coeffs, delta, theta, N, dim, 2)
% end FD check for coefficients


yinterp = zeros(3,res,res);
were_there_failed_processes = 0;

% start point for optimization
q0 = samplocs(:,1); % use first sample loc in first run
                    % for the other runs, use qstar from previous 
                    % run as starting point
error_field = zeros(res,res);

disp("Running interpolation over trial space, please wait.")
tic;
% count number of iterations in gradient descent
count = zeros(res,res);
for j=1:res
    for k=1:res
        wstar = [wspace1_trial(j); wspace2_trial(k)];
        
        % true function value:
        [q_true, d2q_true, d2q_true] = testfun_handle(wspace1_trial(j), wspace2_trial(k));
        
        % interpolate the coefficient funcs
        weights = zeros(N,1);
        for l=1:N
            weights(l) = GEK_interp(wstar, Wlocs, v_GEK_array(:,l), beta_array(l), theta, 1, 2);
        end
        % Do the weight functions sum up to 1? Uncomment the next line.
        %disp(['sum of weights = ', num2str(sum(weights))])

        % find weighted center
        [yinterp(:,j,k), count(j,k), fail] = findCenter_Sphere(samplocs, weights, q0, tau);

        % record error
        error_field(j,k) = norm(yinterp(:,j,k) - q_true, 2);

        % use interpolated value from previous run for next run
        q0 = yinterp(:,j,k);
        were_there_failed_processes = were_there_failed_processes + fail;
    end
end
t_interp = toc;
disp(['data interpolation: ', num2str(t_interp), 's'])




disp([' ']);
disp(['*****************************************']);
disp([' Number of  failures in -findCenter- : ', num2str(were_there_failed_processes)])
disp([' Av. number of iters in -findCenter- : ', num2str(sum(sum(count))/(res*res))])
disp([' Max number of iters in -findCenter- : ', num2str(max(max(count)))])
disp([' Av. interpolation error             : ', num2str(sum(sum(error_field))/(res*res))])
disp([' Max interpolation error             : ', num2str(max(max(error_field)))])
disp(['*****************************************']);


% plot errors:
[Xerrors, Yerrors] = meshgrid(wspace1_trial, wspace2_trial);

subplot(1,2,2)
surf(Xerrors, Yerrors, error_field)
xlabel('$\omega_1$')
ylabel('$\omega_2$')
zlabel('Errors $\|f(\omega) - \hat f(\omega)\|_2$')
colormap gray;
title('Interpolation errors: BHI')
brighten(0.8)

%***************************************
disp([' ']);
disp(['**************************************************']);
disp('Perform FD check for interpolated manifold function')
disp(['**************************************************']);
% reset start point
q0 = mean/norm(mean);
d1_error = zeros(n1*n2,1);
d2_error = zeros(n1*n2,1);
for j=1:n1*n2
    index = j;
    [d1ystar, d2ystar] = FD_check_RiemannInterp(samplocs, Wlocs, q0, tau, v_GEK_array, beta_array, theta, N, dim, delta, index, 'S2');
    %[d1ystar,d1samplocs(:,index), d2ystar, d2samplocs(:,index)]

    d1_error(j) = norm(d1ystar - d1samplocs(:,index))/norm(d1samplocs(:,index));
    d2_error(j) = norm(d2ystar - d2samplocs(:,index))/norm(d1samplocs(:,index));
    disp(['Error( d/d1 f(wj) - d/d1 hat f(wj)) at index j =', num2str(index), ': ',...
        num2str(d1_error(j))])

    disp(['Error( d/d2 f(wj) - d/d2 hat f(wj)) at index j =', num2str(index), ': ',...
        num2str(d2_error(j))])
end
av_FD_err1 = sum(d1_error)/(n1*n2)
av_FD_err2 = sum(d2_error)/(n1*n2)
