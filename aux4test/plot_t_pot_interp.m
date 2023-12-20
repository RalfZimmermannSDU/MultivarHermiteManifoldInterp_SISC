%
% This script plots the tea pot figure under the action of 
% a rotation matrix.
% In the paper
% "MULTIVARIATE HERMITE INTERPOLATION ON RIEMANNIAN MANIFOLDS",
% Zimmermann/Bergmann,
% this correspinds to FIGURE 7.
%
%
% The rotation matrices must be available in form of an
% array of matrices called "refs_mats" in the workspace
%
% To get the list of reference matrices, first run the script
%    "Hermite_RiemannBary_SO.m"
%
% 
%
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