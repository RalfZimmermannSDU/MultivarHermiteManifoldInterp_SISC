function [Wlocs, samplocs, d1samplocs, d2samplocs, pstar, logs_array, d1_logs_array, d2_logs_array] = create_sample_data_SO3_tang(test_fun, Wspace, N, dim)
%
% This function creates the manifold-valued sample data,
% manifold points and tangent vectors.
% !!PROVIDES PREPROCESSING FOR TANGENT SPACE INTERPOLATION!!
% (for preprocessing for barycentric interp, see 
% "create_sample_data_SO3.m")
%
% Function is hard-coded to work for data on SO(3)
%
% INPUTS
% test_fun : function handle for test function, 
%            must provide samples plus derivative
% wspacej  : sample locations along j-th coordinate direction. j=1,2
% N        : number of sample points
% dim      : dimension (actually hard coded dim = 2)
%
% OUTPUTS
% Wlocs      : sample locations in parameter space
% samplocs   : sample points on manifold
% djsamplocs : partial manifold derivative in coordinate direction j=1,2
%              at sample locations
% djcoeffs   : j-partial derivatives of weight functions, j=1,2
%
%
%

% dimension of SO-matrices
d = 3;

Wlocs      = zeros(dim,N);
samplocs   = zeros(d,d,N);
d1samplocs = zeros(d,d,N);
d2samplocs = zeros(d,d,N);

Wlocs = Wspace;

% compute sample locs and derivatives
for j=1:N
        % sample locs + derivatives
        [qw, d1qw, d2qw]  = test_fun(Wspace(1,j), Wspace(2,j));
        % store location in parameter space
        samplocs(:,:,j)   = qw;
        d1samplocs(:,:,j) = d1qw;
        d2samplocs(:,:,j) = d2qw;
end

%-----------------------------------------
% compute tangent space base point
%-----------------------------------------
weights = (1/N)*ones(N,1);      % equal weights
q0 = samplocs(:,:,floor(N/2));  % initial guess
tau = 1.0e-8;
[pstar, count, fail] =  findCenter_SOn(samplocs, weights, q0, tau);

%-----------------------------------------
% compute tangent space images
%-----------------------------------------
max_dist =0.0;

d1coeffs   = zeros(N,N);
d2coeffs   = zeros(N,N);
logs_array = zeros(d,d,N);   % storage space for Riemann logs
    
%map all data to T_(pstar)S
for i = 1:N     
    Logi = log_SOn(pstar,samplocs(:,:,i));
    dist_ij = norm(Logi,"fro");
    if dist_ij > max_dist
        max_dist = dist_ij;
    end
    logs_array(:,:,i) = Logi;
end

disp(['Max dist in log calcs: ', num2str(max_dist/pi), 'pi'])

% compute the tangent space images of the sampled derivatives
% hard coded for 2D

d1_logs_array = zeros(d,d,N);
d2_logs_array = zeros(d,d,N);

for i=1:N
    % data for central difference approximation
    h = 1.0e-5;
    Tplus = exp_SOn(samplocs(:,:,i), h*d1samplocs(:,:,i));
    fplus = log_SOn(pstar, Tplus);
    Tminus= exp_SOn(samplocs(:,:,i), (-h)*d1samplocs(:,:,i));
    fminus= log_SOn(pstar, Tminus);
    % central difference approximation   
    d1_logs_array(:,:,i) = (1.0/(2*h))*(fplus - fminus);

    % for second dim
    Tplus = exp_SOn(samplocs(:,:,i), h*d2samplocs(:,:,i));
    fplus = log_SOn(pstar, Tplus);
    Tminus= exp_SOn(samplocs(:,:,i), (-h)*d2samplocs(:,:,i));
    fminus= log_SOn(pstar, Tminus);
    % central difference approximation   
    d2_logs_array(:,:,i) = (1.0/(2*h))*(fplus - fminus);
end


return;
end