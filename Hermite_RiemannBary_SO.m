%
%
% This script performs the experiment from Section 5.2 of the paper
% "MULTIVARIATE HERMITE INTERPOLATION ON RIEMANNIAN MANIFOLDS",
%  Zimmermann/Bergmann
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
n1    = 9  % sample resolution in coordinate direction 1
n2    = 7  % sample resolution in coordinate direction 1
res   = 76 % resolution for trial space
delta = 1.0e-5     % delta for finite difference tests
tau   = 1.0e-6     % convergence threshold for Riemannian optimization
index = 1          % check interpolation condition for sample point "index"
theta = [0.5,0.5]; % correlation weights
%
%
do_midpoint_trials = 1; % boolean: if set to "1", resolution an
                        % trial space is adjusted so that 
                        % the tea pot figure on six midpoint trial points
                        % can be created with the script
                        % aux4test/plot_t_pot_interp
                        % hard coded for n1=n2, res = n1-1
if do_midpoint_trials
    n2  = n1;
    res = n1-1;
end
%
%
% create sample points 
%
%
%
a = -.5;
b = .5;
%

sample_meth = input("Sampling method: 1(full fac),2(cheby),3(rLHC): ");

if sample_meth ==1
    wspace1 = linspace(a,b,n1); %(b-a).*rand(n1,1) + a; 
    wspace2 = linspace(a,b,n2); %(b-a).*rand(n2,1) + a;
    % create an initial sample locs array of the form
    % w11, w21, w31, ..., wn1
    % w12, w22, w32, ..., wn2
    [W1, W2] = meshgrid(wspace1, wspace1);
    Wspace = [W1(:)'; W2(:)'];
elseif sample_meth ==2
    wspace1 = ChebyRoots(a, b, n1);
    wspace2 = ChebyRoots(a, b, n2);
    [W1, W2] = meshgrid(wspace1, wspace1);
    Wspace = [W1(:)'; W2(:)'];
elseif sample_meth ==3
    % random latin hypercube
    X = rlh(n1*n2, 2);
    Wspace = (b-a).*X' + a;
end



% Execute the following code block 
% for doing midpoint trials and creating 
% the reference data set for producing the tea pot figure.
if do_midpoint_trials
    % this code bock creates trial points 
    % exactly in the midpoints of the two-D sample grid
    for j=1:length(wspace1)-1
        wspace1_trial(j) = wspace1(j) + 0.5*(wspace1(j+1) -wspace1(j));
    end
    %
    for j=1:length(wspace2)-1
        wspace2_trial(j) = wspace2(j) + 0.5*(wspace2(j+1) -wspace2(j));
    end
else
    % standard trial space
    wspace1_trial = linspace(a,b,res);
    wspace2_trial = linspace(a,b,res);
end

dim = 2; % 2D sample space

N = n1*n2;

%testfun_handle = @(x,y) testfun1_S2(x,y);
testfun_handle = @(x,y) testfun_SO3(x,y);


tic;
% create all data required for Hermite interpolation
% choose test function
[Wlocs, samplocs, d1samplocs, d2samplocs, d1coeffs, d2coeffs] = create_sample_data_SO3(testfun_handle, Wspace, N, dim);
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
    [v_GEK, beta]    = setup_GEK_interp(Wlocs, Yphi(:,k), theta);
    v_GEK_array(:,k) = v_GEK;
    beta_array(k)    = beta;
end
t_preprocessing = toc;
disp(['data preprocessing: ', num2str(t_preprocessing), 's'])

% dimension of matrices
d = size(samplocs, 1);

% compute some q as a candidate for the barycenter
mean = zeros(d,d);
for l=1:N 
    Delta_l = log_SOn(samplocs(:,:,floor(N/2)), samplocs(:,:,l));
    mean = mean + Delta_l;
end
mean = mean/N;
qmean = exp_SOn(samplocs(:,:,floor(N/2)), mean);


% FD check for coefficients
FD_test_coeffs(Wlocs, d1coeffs, d2coeffs, delta, theta, N, dim, 2)
% end FD check for coefficients



yinterp  = zeros(d,d,res,res);
ref_mats = zeros(d,d,res,res);
were_there_failed_processes = 0;

% start point for optimization
q0 = samplocs(:,:,1); % use first sample loc in first run
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
        ref_mats(:,:,j,k) = q_true;

        % interpolate the coefficient funcs
        weights = zeros(N,1);
        for l=1:N
            weights(l) = GEK_interp(wstar, Wlocs, v_GEK_array(:,l), beta_array(l), theta, 1, 2);
        end
        % Do the weight functions sum up to 1?
        %disp(['sum of weights = ', num2str(sum(weights))])

        % find weighted center
        [yinterp(:,:,j,k), count(j,k), fail] = findCenter_SOn(samplocs, weights, q0, tau);

        % record error
        % Fro-norm of an SO3 matrix is sqrt(3)
        error_field(j,k) = norm(yinterp(:,:,j,k) - q_true, 'fro')/sqrt(3);

        % use interpolated value from previous run for next run
        q0 = yinterp(:,:,j,k);
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

figure
surf(Xerrors, Yerrors, error_field)
xlabel('$\omega_1$','interpreter','latex')
ylabel('$\omega_2$','interpreter','latex')
zlabel('Errors $\|f(\omega) - \hat f(\omega)\|_2$','interpreter','latex')
colormap gray;
title('Interpolation errors: BHI', 'interpreter','latex')
brighten(0.8)
%
%
%
%


if 1
    %***************************************
    disp([' ']);
    disp(['**************************************************']);
    disp('Perform FD check for interpolated manifold function')
    disp(['**************************************************']);
    % reset start point
    %q0 = qmean;
    q0 = samplocs(:,:,1);

    d1_error = zeros(n1*n2,1);
    d2_error = zeros(n1*n2,1);
    for j=1:n1*n2
        index = j;
        [d1ystar, d2ystar] = FD_check_RiemannInterp(samplocs, Wlocs, q0, tau, v_GEK_array, beta_array, theta, N, dim, delta, index, 'SO');
        %[d1ystar,d1samplocs(:,index), d2ystar, d2samplocs(:,index)]

        d1_error(j) = norm(d1ystar - d1samplocs(:,:,index), 'fro')/norm(d1samplocs(:,:,index), 'fro');
        d2_error(j) = norm(d2ystar - d2samplocs(:,:,index), 'fro')/norm(d1samplocs(:,:,index), 'fro');
        disp(['Error( d/d1 f(wj) - d/d1 hat f(wj)) at index j =', num2str(index), ': ',...
            num2str(d1_error(j))])

        disp(['Error( d/d2 f(wj) - d/d2 hat f(wj)) at index j =', num2str(index), ': ',...
            num2str(d2_error(j))])
    end
    av_FD_err1 = sum(d1_error)/(n1*n2)
    av_FD_err2 = sum(d2_error)/(n1*n2)
end