function [dist] = dist_SOn(X,Y)
% Riemannian distance on SOn
dist = norm(logm(X'*Y), 'fro');
end