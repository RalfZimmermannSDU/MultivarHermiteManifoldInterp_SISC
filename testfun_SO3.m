function [qw, d1qw, d2qw] = testfun_SO3(w1,w2)
% 
% synthetic test function on SO3
%
a = 4*pi;
%
f = [w1^2 + .5*w2; sin(a*(w1^2 + w2^2)); w1 + w2^2];


X = [[0, f(1), f(2)]; [-f(1), 0,  f(3)]; [-f(2), -f(3), 0]];
d1X = [[0, 2*w1, cos(a*(w1^2 + w2^2))*a*2*w1]; [-2*w1, 0,  1]; [-cos(a*(w1^2 + w2^2))*a*2*w1, -1, 0]];
d2X = [[0, .5, cos(a*(w1^2 + w2^2))*a*2*w2]; [-.5, 0,  2*w2]; [-cos(a*(w1^2 + w2^2))*a*2*w2, -2*w2, 0]];

% compute fun values and derivatives with Mathias' theorem
T1 = [[X, d1X]; [zeros(size(X)), X]];
T2 = [[X, d2X]; [zeros(size(X)), X]];

expT1 = expm(T1);
expT2 = expm(T2);

qw = expT1(1:3,1:3);
%partial derivatives

d1qw = expT1(1:3,4:6);
d2qw = expT2(1:3,4:6);

% do we have proper sample points on SO3
check2 = qw'*qw -eye(3);
% are the derivatives proper tangent vectors?
check1 = [d1qw'*qw + qw'*d1qw, d2qw'*qw + qw'*d2qw];
if norm(check1)> 1.0e-13 || norm(check2)> 1.0e-13
    disp('ARGH  (in testfun_SO3)');
    return;
end
%
return
%
end