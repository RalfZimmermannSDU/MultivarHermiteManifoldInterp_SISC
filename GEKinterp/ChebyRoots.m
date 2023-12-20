function [chebyroots] = ChebyRoots(a, b, n)
%--------------------------------------------------------------------------
% compute n chebychev roots in the interval [a,b]
%--------------------------------------------------------------------------

chebyroots = zeros(n,1);
    
% in [-1,1]-interval
for i=1:n
    % filling the array from the rear gives the roots in ascending order
    chebyroots(n+1 -i) = cos((2*i-1)*pi/(2*n));
end
% translate to [a,b]
chebyroots = 0.5*(b-a)*chebyroots + 0.5*(b+a);
return 
%END ----------------------------------------------------------------------
end
