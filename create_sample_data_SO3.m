function [Wlocs, samplocs, d1samplocs, d2samplocs, d1coeffs, d2coeffs] = create_sample_data_SO3(test_fun, Wspace, N, dim)
%
% This function creates the manifold-valued sample data,
% manifold points and tangent vectors.
% !!PROVIDES PREPROCESSING FOR BARYCENTRIC INTERPOLATION!!
% (for preprocessing for TANGENT SPACE interp, see 
% "create_sample_data_SO3_tang.m")
%
% Function is hard-coded to work for data on SO(3).
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

Wlocs = zeros(dim,N);
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
% compute coeffs of weight functions phi 
%-----------------------------------------
max_dist = 0.0;
%for k=1:dim  % here manual for 2D
    d1coeffs = zeros(N,N);
    d2coeffs = zeros(N,N);
    for j=1:N
        logs_array_vec = zeros(d*d,N);   % storage space for Riemann logs
        %map all data to TpjS, store vectorized matrices.
        for i = 1:N     
            if i~=j
                Logi = log_SOn(samplocs(:,:,j),samplocs(:,:,i));
                dist_ij = norm(Logi,"fro");
                if dist_ij > max_dist
                    max_dist = dist_ij;
                end
                logs_array_vec(:,i) = Logi(:);
            end
        end
        % the tangent space of S2 is of dim 2. 
        % Do the logs span the full tangent space?
        
        % The current sample loc is pj  =   samplocs(:,:,j)
        % The derivatives are       v1j = d1samplocs(:,:,j)
        %                           v2j = d2samplocs(:,:,j)

        % to compute the Hessian is not necessary, since it is the identity
        % in the only cases that are needed.
        % Hess_vkj = Hess_dist(p1,vkj,p1) = vkj

        % Hence the task is to represent vkj as a linear combination of
        % the logs

        % -->solve under-determined(?) least-squares system
        % remove zero column
        X = [logs_array_vec(:,1:j-1), logs_array_vec(:,j+1:N)];

        
        % -->solve under-determined(?) least-squares system

        [VT, Sig, UT] = svd(X',0);
        % perform calculations in U-coordinates
        Sig_vec = diag(Sig);
        rankL = sum(Sig_vec > 1.0e-10);
        % now VT*sig*UT' = X' <=> UT*Sig*VT' = X.
        %remove zero sing vals
        Sig_inv = diag(1.0./Sig_vec(1:rankL));
        %Ur = UT(1:rankL,:)';
        Ur = UT(:,1:rankL);

        % Ur  is a basis for the column space of X
        % need to solve Ur'*X*c = Ur'*v
        Xr = [Ur'*X; ones(1,size(X,2))];
        v1jr = [Ur'*reshape(d1samplocs(:,:,j), d*d,1); 0];
        v2jr = [Ur'*reshape(d2samplocs(:,:,j), d*d,1); 0];
        dcoeffs_j = linsolve(Xr*Xr', [v1jr, v2jr]);
        dcoeffs_j = Xr'*dcoeffs_j;
        d1coeffs_j = dcoeffs_j(:,1);
        d2coeffs_j = dcoeffs_j(:,2);
        %% does c solve Xc = vkj?
        %check_solutiond1 = norm(reshape(d1samplocs(:,:,j), d*d,1) - X*d1coeffs_j)/norm(reshape(d1samplocs(:,:,j), d*d,1))
        %check_solutiond2 = norm(reshape(d2samplocs(:,:,j), d*d,1) - X*d2coeffs_j)/norm(reshape(d2samplocs(:,:,j), d*d,1))
        %
        %% alternative solution via SVD and Schur complement inversion
        %km1 = size(X,2)
        %size(VT)
        %Vr = VT(:,1:rankL);
        %Ur = UT(:,1:rankL);
        
        %Sr = diag(Sig_vec(1:rankL))
        %norm(X - Ur*Sr*Vr')
        %v1jr = Ur'*reshape(d1samplocs(:,:,j), d*d,1);
        %xtmp=Vr'*ones(size(X,2), 1);
        %vtmp = Vr*((1.0./Sig_vec(1:rankL)).*v1jr);
        %S = km1 - xtmp'*xtmp
        %(1/S)
        %c1 = (1/S)*( S*vtmp + (ones(1,km1)*vtmp)*(Vr*Vr'-eye(km1))*ones(km1,1))
        %% End alternative------------------


        % copy to coefficients array, skip diagonal entry
        d1coeffs(1:(j-1),j) = d1coeffs_j(1:j-1);
        d1coeffs(j+1:N,j)   = d1coeffs_j(j:N-1);
        d2coeffs(1:(j-1),j) = d2coeffs_j(1:j-1);
        d2coeffs(j+1:N,j)   = d2coeffs_j(j:N-1);
    end
%end

disp(['Max dist in log calcs: ', num2str(max_dist/pi), 'pi'])

return;
end