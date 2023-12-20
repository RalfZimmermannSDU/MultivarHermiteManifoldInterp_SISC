function [logS] = log_sphere(q,p)
% Riemann log on sphere
%
% if q and p are numerically identical, return 0
if norm(q-p,2) < 100*eps;
    logS = zeros(size(q));
    return;
end
%
%
inner_q_p = q'*p;
% inner product must be of abs < 1; => filter numerical errors
if abs(inner_q_p)>1.0
    inner_q_p = sign(inner_q_p)*1.0;
end

logS = (acos(inner_q_p)/norm(p-(inner_q_p)*q,2))*(p-(inner_q_p)*q);

return;
end