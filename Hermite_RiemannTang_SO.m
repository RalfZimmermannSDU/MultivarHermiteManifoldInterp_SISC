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
n1    = 7  % sample resolution in coordinate direction 1
n2    = 7  % sample resolution in coordinate direction 1
res   = 76 % resolution for trial space
delta = 1.0e-5     % delta for finite difference tests
tau   = 1.0e-8     % convergence threshold for Riemannian optimization
index = 1          % check interpolation condition for sample point "index"
theta = [0.5,0.5]; % correlation weights
do_checks = 1;     % Boolean variable for doing extra checks

%
%
% create sample points 
%
%
%
a =-.5;
b = .5;
%

sample_meth = input("Sampling method = 1 (full fac),2(cheby),3(rLHC): ");

if sample_meth ==1
    % full factorial sampling in 2D
    wspace1 = linspace(a,b,n1);  
    wspace2 = linspace(a,b,n2); 
    % create an initial sample locs array of the form
    % w11, w21, w31, ..., wn1
    % w12, w22, w32, ..., wn2
    [W1, W2] = meshgrid(wspace1, wspace1);
    Wspace = [W1(:)'; W2(:)'];
elseif sample_meth ==2
    % 2D sampling on Chebychev roots
    wspace1 = ChebyRoots(a, b, n1);
    wspace2 = ChebyRoots(a, b, n2);
    [W1, W2] = meshgrid(wspace1, wspace1);
    Wspace = [W1(:)'; W2(:)'];
elseif sample_meth ==3
    % random latin hypercube
    X = rlh(n1*n2, 2);
    % map to data range
    Wspace = (b-a).*X' + a;
end


% trial space: these are the locations, 
% where the interpolation method is assessed
wspace1_trial = linspace(a,b,res);
wspace2_trial = linspace(a,b,res);


% the next for loop creates trial samples at the "center cells" along a
% diagonal in a 2D sample space
% This is used for the tee-pot experiment in the manuscript
% "Multivariate Hermite Interpolation on Manifolds
test_along_diag = 0;
if test_along_diag
    for j=1:length(wspace1)-1
        wspace1_trial(j) = wspace1(j) + 0.5*(wspace1(j+1) -wspace1(j));
    end
    for j=1:length(wspace2)-1
        wspace2_trial(j) = wspace2(j) + 0.5*(wspace2(j+1) -wspace2(j));
    end
end

dim = 2;   % 2D sample space, method is hard-coded for this case.

N = n1*n2; % total number of samples

% specify test function
testfun_handle = @(x,y) testfun_SO3(x,y);

% Preprocessing
% create all data required for Hermite interpolation
tic;
[Wlocs, samplocs, d1samplocs, d2samplocs, pstar, logs_array, d1_logs_array, d2_logs_array] = create_sample_data_SO3_tang(testfun_handle, Wspace, N, dim);
% dimension of matrices
d = size(samplocs, 1);


% set up sample data vecs for interpolation
Yphi = eye(N*(dim+1),N*(dim+1));
% each column j of Yphi contains
%   * the sample values of the coeff funcs
%   * the partial derivatives in direction d1 of the coeff funcs
%   * the partial derivatives in direction d2 of the coeff funcs
% for the linear combination
% The interpolant is of the form
% yhat = \sum_j \lambda_j(w) logs_array(j) + \sum_ij \mu_ij(w) di_logs_array(j)
%
% the sample data vector is one column of Yphi, namely
% Yphi(:, k) = (lambda_k(w1),...,lambda_k(wn),
% d1lambda_k(w1),...,d1lambda_k(wn), ...)

% Next step: Interpolate the coefficient functions 

% set up Kriging predictor for each Yphi(:,k)
v_GEK_array = zeros(N*(dim+1), N*(dim+1));
beta_array  = zeros(N*(dim+1),1);
for k=1:N*(dim+1)
    [v_GEK, beta] = setup_GEK_interp(Wlocs, Yphi(:,k), theta);
    % outputs are:
    %  * the GEK vector v_GEK
    %  * the regression coefficient beta
    v_GEK_array(:,k) = v_GEK;
    beta_array(k)    = beta;
end
t_preprocessing = toc;
disp(['data preprocessing: ', num2str(t_preprocessing), 's'])

% storage array for interpolation results
yinterp = zeros(d,d,res,res);
% storage array for interpolation errors
error_field = zeros(res,res);

%
disp("Running interpolation over trial space, please wait.")
tic;
for j=1:res
    for k=1:res
        wstar = [wspace1_trial(j); wspace2_trial(k)];
        
        % compute true function value for comparison purposes
        [q_true, d2q_true, d2q_true] = testfun_handle(wspace1_trial(j), wspace2_trial(k));
        
        % interpolate the coefficient funcs
        weights = zeros(N*(dim+1),1);
        for l=1:N*(dim+1)
            weights(l) = GEK_interp(wstar, Wlocs, v_GEK_array(:,l), beta_array(l), theta, 1, 2);
        end

        if 0
            % Do the weight functions sum up to 1?
            disp(['sum of sample weights = ', num2str(sum(weights(1:N)))])
        end

        % interpolate the tangent space  matrices
        y_tang = zeros(d,d);
        for l=1:N
            y_tang = y_tang + weights(l)*logs_array(:,:,l) + weights(N+l)*d1_logs_array(:,:,l) + ...
                     weights(2*N+l)*d2_logs_array(:,:,l);
        end
        % map y_tang to manifold
        yinterp(:,:,j,k) = exp_SOn(pstar, y_tang);

        % record error
        % Fro-norm of an SO3 matrix is sqrt(3)
        error_field(j,k) = norm(yinterp(:,:,j,k) - q_true, 'fro')/sqrt(3);
    end
end
t_interp = toc;
disp(['data interpolation tangent space: ', num2str(t_interp), 's'])

disp([' ']);
disp(['*****************************************']);
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
title('Interpolation errors: TSHI', 'interpreter','latex')
brighten(0.8)

%
%
%
%



if do_checks
    %***************************************
    disp([' ']);
    disp(['**************************************************']);
    disp('Perform FD check for interpolated manifold function')
    disp(['**************************************************']);
    
    d1_error = zeros(n1*n2,1);
    d2_error = zeros(n1*n2,1);
    for j=1:n1*n2
        index = j;
    
        wstar    = Wlocs(:,index);
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
        y_tang_1p = zeros(d,d);
        y_tang_1m = zeros(d,d);
        y_tang_2p = zeros(d,d);
        y_tang_2m = zeros(d,d);
        for l=1:N
            y_tang_1p = y_tang_1p + weights_1p(l)*logs_array(:,:,l) + weights_1p(N+l)*d1_logs_array(:,:,l) + ...
                        weights_1p(2*N+l)*d2_logs_array(:,:,l);
            y_tang_1m = y_tang_1m + weights_1m(l)*logs_array(:,:,l) + weights_1m(N+l)*d1_logs_array(:,:,l) + ...
                        weights_1m(2*N+l)*d2_logs_array(:,:,l); 
            y_tang_2p = y_tang_2p + weights_2p(l)*logs_array(:,:,l) + weights_2p(N+l)*d1_logs_array(:,:,l) + ...
                        weights_2p(2*N+l)*d2_logs_array(:,:,l);
            y_tang_2m = y_tang_2m + weights_2m(l)*logs_array(:,:,l) + weights_2m(N+l)*d1_logs_array(:,:,l) + ...
                        weights_2m(2*N+l)*d2_logs_array(:,:,l);
        end
    
        % map to manifold
        y_1p = exp_SOn(pstar, y_tang_1p);
        y_1m = exp_SOn(pstar, y_tang_1m);
        y_2p = exp_SOn(pstar, y_tang_2p);
        y_2m = exp_SOn(pstar, y_tang_2m);
    
        d1ystar = (1.0/(2*delta)) *(y_1p - y_1m);
        d2ystar = (1.0/(2*delta)) *(y_2p - y_2m);
    
        d1_error(j) = norm(d1ystar - d1samplocs(:,:,index), 'fro')/norm(d1samplocs(:,:,index), 'fro');
        d2_error(j) = norm(d2ystar - d2samplocs(:,:,index), 'fro')/norm(d1samplocs(:,:,index), 'fro');
        %disp(['Error( d/d1 f(wj) - d/d1 hat f(wj)) at index j =', num2str(index), ': ',...
        %    num2str(d1_error(j))])
        
        %disp(['Error( d/d2 f(wj) - d/d2 hat f(wj)) at index j =', num2str(index), ': ',...
        %    num2str(d2_error(j))])
    end
    av_FD_err1 = sum(d1_error)/(n1*n2)
    av_FD_err2 = sum(d2_error)/(n1*n2)
end
