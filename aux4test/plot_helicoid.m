% plot Helicoid
%
% This function opens a subplot figure with Helicoid on the left
%
%
function [] = plot_helicoid(uspace, vspace)

% user input

do_sphere_plot = input('Do you want to plot the samples on the unit sphere? 0/1 :' )

% component functions of Helicoid
funx = @(u,v)  sinh(u).*sin(v);
funy = @(u,v) -sinh(u).*cos(v);
funz = @(u,v)  -v;

a =-0.25*pi;
b = 0.25*pi;
ploteps = 0.5; % extend plot range beyond samples

aplot = a-ploteps;
bplot = b+ploteps;

[ugrid, vgrid] = meshgrid(uspace, vspace);

% base coords
X = [funx(ugrid,vgrid)];
Y = [funy(ugrid,vgrid)];
Z = [funz(ugrid,vgrid)];

% normal vec
U = [(1.0./(exp(2.*ugrid) +1)).*(2.*exp(ugrid).*cos(vgrid))];
V = [(1.0./(exp(2.*ugrid) +1)).*(2.*exp(ugrid).*sin(vgrid))];
W = [(exp(2*ugrid) - 1)./(exp(2.*ugrid) +1)];

% partial derivatives
d1factor1 = (-2.*exp(2.*ugrid)./((exp(2.*ugrid)+1).^2));
d1factor2 = (1.0./(exp(2.*ugrid)+1));

d1U = d1factor1.*(2.*exp(ugrid).*cos(vgrid)) + ...
      d1factor2.*(2.*exp(ugrid).*cos(vgrid));
d1V = d1factor1.*(2.*exp(ugrid).*sin(vgrid)) + ...
      d1factor2.*(2.*exp(ugrid).*sin(vgrid));
d1W = d1factor1.*(exp(2.*ugrid) - 1) + ...
      d1factor2.*(2.*exp(2.*ugrid));

d2factor = (1.0./(exp(2.*ugrid)+1));
d2U = d2factor.*(-2.*exp(ugrid).*sin(vgrid));
d2V = d2factor.*(2.*exp(ugrid).*cos(vgrid));
d2W = d2factor.*0.0;




% additional plot of samples on unit sphere
if do_sphere_plot
    figure;
    % plot sphere
    sphere;
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    alpha(0.3);
    title('Sample points on sphere')
    hold on 
    X0 = zeros(size(U));
    quiver3(X0,X0,X0,U,V,W, 1, 'k', 'AutoScale','off', 'LineWidth', 2)
end

figure;
% this is the Helicoid lus sample data
subplot(1,2,1);
fsurf(funx,funy,funz,[ aplot bplot aplot bplot]) 
camlight
colormap(gray)
alpha(0.5)
title('Sample data: normal vectors on helicoid')
hold on


quiver3(X,Y,Z,U,V,W, 1, 'k', 'AutoScale','off', 'LineWidth', 3)
hold on
quiver3(X,Y,Z,d1U,d1V,d1W, 1, 'r', 'AutoScale','off', 'LineWidth', 1)
hold on
quiver3(X,Y,Z,d2U,d2V,d2W, 1, 'b', 'AutoScale','off', 'LineWidth', 1)
hold on
end
