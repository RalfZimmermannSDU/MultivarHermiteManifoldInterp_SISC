function [LogXY] = log_SOn(X, Y)
%
% Riemannian log on SOn
%
% returns Log_X(Y)

LogXY = X*logm(X'*Y);

%
end