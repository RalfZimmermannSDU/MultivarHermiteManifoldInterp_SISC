function [ExpXD] = exp_SOn(X, D)
%
% Riemannian exp on SOn
%
% returns Exp_X(D)

% for numerical safety: make input skew
XTD = X'*D;
XTD = 0.5*(XTD-XTD');

ExpXD = X*expm(XTD);

return;
%
end