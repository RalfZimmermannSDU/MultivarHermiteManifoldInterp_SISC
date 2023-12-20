function [Wlocs, samplocs, d1samplocs, d2samplocs, d1coeffs, d2coeffs] = create_sample_data_Sphere(test_fun, wspace1, wspace2, N, dim)
% This function creates the manifold-valued sample data,
% i.e. manifold points and tangent vectors.
% !!PROVIDES PREPROCESSING FOR BARYCENTRIC INTERPOLATION!!
% (for preprocessing for tangent space interp, see 
% "create_sample_data_Sphere_tang.m"
%
% Function is hard-coded to work for data on the unit sphere
%
% INPUTS
% test_fun : function handle for test function from para space to sphere, 
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
Wlocs = zeros(dim,N);
samplocs   = zeros(3,N);
d1samplocs = zeros(3,N);
d2samplocs = zeros(3,N);

n1 = length(wspace1);
n2 = length(wspace2);


% compute sample locs and derivatives
for j=1:n1
    for k=1:n2
        % sample locs + derivatives
        [qw, d1qw, d2qw] = test_fun(wspace1(j), wspace2(k));
        % store location in parameter space
        Wlocs(:,k+n2*(j-1)) = [wspace1(j); wspace2(k)];
        samplocs(:,k+n2*(j-1)) = qw;
        d1samplocs(:,k+n2*(j-1)) = d1qw;
        d2samplocs(:,k+n2*(j-1)) = d2qw;
    end
end


%-----------------------------------------
% compute coeffs of weight functions phi 
%-----------------------------------------

%for k=1:dim  % here manual for 2D
    d1coeffs = zeros(N,N);
    d2coeffs = zeros(N,N);
    for j=1:N
        logs_array = zeros(3,N);   % storage space for logs
        %map all data to TpjS
        for i = 1:N     
            if i~=j
                logs_array(:,i) = log_sphere(samplocs(:,j),samplocs(:,i));
            end
        end
        % the tangent space of S2 is of dim 2. 
        % Do the logs span the full tangent space?
        
        % The current sample loc is pj  =   samplocs(:,j)
        % The derivatives are       v1j = d1samplocs(:,j)
        %                           v2j = d2samplocs(:,j)

        % to compute the Hessian is not necessary, since it is the identity
        % in the only cases that are needed.
        % Hess_vkj = Hess_dist(p1,vkj,p1) = vkj

        % Hence, the task is to represent vkj as a linear combination of
        % the logs

        %remove zero column
        X = [logs_array(:,1:j-1), logs_array(:,j+1:N)];

        
        % -->solve under-determined(?) least-squares system
        [VT, Sig, UT] = svd(X',0);
        % perform calculations in U-coordinates
        Sig_vec = diag(Sig);
        rankL = sum(Sig_vec > 1.0e-10);
        % now VT*sig*UT' = X' <=> UT*Sig*VT' = X.
        %remove zero sing vals
        Sig_inv = diag(1.0./Sig_vec(1:rankL));
        Ur = UT(:,1:rankL);
        % Ur  is a basis for the column space of X
        % need to solve Ur'*X*c = Ur'*v
        Xr = [Ur'*X; ones(1,size(X,2))];
        v1jr = [Ur'*d1samplocs(:,j); 0];
        v2jr = [Ur'*d2samplocs(:,j); 0];
        dcoeffs_j = linsolve(Xr*Xr', [v1jr, v2jr]);
        dcoeffs_j = Xr'*dcoeffs_j;
        d1coeffs_j = dcoeffs_j(:,1);
        d2coeffs_j = dcoeffs_j(:,2);
        %% does c solve Xc = vkj?
        %check_solutiond1 = norm(d1samplocs(:,j) - X*d1coeffs_j)/norm(d1samplocs(:,j))
        %check_solutiond2 = norm(d2samplocs(:,j) - X*d2coeffs_j)/norm(d2samplocs(:,j))
        

        % copy to coefficients array, skip diagonal entry
        d1coeffs(1:(j-1),j) = d1coeffs_j(1:j-1);
        d1coeffs(j+1:N,j)   = d1coeffs_j(j:N-1);
        d2coeffs(1:(j-1),j) = d2coeffs_j(1:j-1);
        d2coeffs(j+1:N,j)   = d2coeffs_j(j:N-1);
    end
%end
return;
end