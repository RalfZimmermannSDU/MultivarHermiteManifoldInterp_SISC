function [expS] = exp_sphere(q,t,v)
% Riemann exponential on sphere
expS = cos(t*norm(v,2))*q + (sin(t*norm(v,2))/norm(v,2))*v;
%
% norm should be one
% add double calculation for numerical robustness
expS = expS/norm(expS,2);
%
return;
%
end