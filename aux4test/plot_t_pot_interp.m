%
% This script plots the tea pot figure under the action of 
% a rotation matrix.
%
% In the paper
% "MULTIVARIATE HERMITE INTERPOLATION ON RIEMANNIAN MANIFOLDS",
% Zimmermann/Bergmann,
% this corresponds to FIGURE 7.
%
% The script is hard coded to illustrate the first 6 matrices along the
% diagonal in the 2D trial space.
% For generating the figure, the rotation matrices must be 
% available in form of an array of matrices called 
% "ref_mats" in the Matlab workspace.
%
% To get the list of reference matrices from the figure,
% first open the script
%    "Hermite_RiemannBary_SO.m",
% set "n1 = 7", 
% set the boolean variable "do_midpoint_trials = 1"
% and run the script.
% 
% Afterwards, run this script.
%

%number of pots to plot
nc = 6; % six reference figures compared with interpolated
        % => 12-tea pots, arrange on (3 x 4)-grid of figures


[verts0, faces, cindex] = teapotGeometry;

figure;
tiledlayout(3,4, 'Padding', 'none', 'TileSpacing', 'compact'); 

count  = 1;
for c = 1:2:2*nc
    subplot(3,4,c);
    
    verts = verts0*ref_mats(:,:,count,count)

    p = patch('Faces',faces,'Vertices',verts,'FaceVertexCData',cindex,'FaceColor','interp', 'FaceAlpha', 1.0);

    p.LineStyle = 'none';      % remove the lines

    view(50,60)   % change the orientation
    axis equal    % make the axes equal and invisible
    %axis([ -3.5 3.5 -3.5 3.5 -0.5 3.5])
    ax = gca;
    ax.XTickLabel = {};
    ax.YTickLabel = {};
    ax.ZTickLabel = {};
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.ZGrid = 'on';
    colormap(bone)
    
    subplot(3,4,c+1);

    verts_rot = verts0*yinterp(:,:,count,count);

    prot = patch('Faces',faces,'Vertices',verts_rot,'FaceVertexCData',cindex,'FaceColor','b', 'FaceAlpha', 0.35)
    prot.LineStyle = 'none';      % remove the lines
    view(50,60)     % change the orientation
    axis equal      % make the axes equal and invisible
    ax = gca;
    ax.XTickLabel = {};
    ax.YTickLabel = {};
    ax.ZTickLabel = {};
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.ZGrid = 'on';
    %axis([ -3.5 3.5 -3.5 3.5 -0.5 3.5])
    count = count +1
end
