function [qw, d1qw, d2qw] = testfun_gauss_S2(w1,w2)
% test function on S^2 c R^3
% Gauss map of Helicoid, stereographic projection of e^(w1+iw2)
%
%
qw = zeros(3,1);
d1qw = zeros(3,1);
d2qw = zeros(3,1);



qw(1) = 2*exp(w1)*cos(w2);
qw(2) = 2*exp(w1)*sin(w2);
qw(3) = exp(2*w1) - 1;

qw = qw/(exp(2*w1) +1);



d1qw = (-2*exp(2*w1)/( (exp(2*w1)+1)^2)) * [2*exp(w1)*cos(w2);2*exp(w1)*sin(w2);exp(2*w1) - 1] ...
       + (1.0/(exp(2*w1)+1)) * [2*exp(w1)*cos(w2);2*exp(w1)*sin(w2); 2*exp(2*w1)];
d2qw =  (1.0/(exp(2*w1)+1))* [-2*exp(w1)*sin(w2);2*exp(w1)*cos(w2); 0.0];

% are the derivatives proper tangent vectors?
check1 = [d1qw'*qw, d2qw'*qw];
if norm(check1)> 1.0e-14
    disp('ARGH  (in testfun_gauss_S2)');
    return;
end
%
%
% for numerical robustness: double-normalize qw
qw = qw/norm(qw);

%
end