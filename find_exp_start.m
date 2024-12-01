function [exponent,int] = find_exp_start(A, f, tol, start)
% Finds the exponential m, which we can use to truncate the Taylor series of the determinant
% A = cell array containing the matrix coefficients
% f = function handling the analytic functions
% tol = tolerance we ask in the approximation of the determinant

if (nargin < 4)
    start = 5;
end

expon = start:150;
s = length(expon);

for k = 1:s
    j = expon(k);
    
    [approx, ~] = trap_rule_fixpow(A, f, j);
    
    if (approx <= tol)
        
         exponent = j;
         int = approx;
         break
    else
        
        exponent = j;
        int = approx;
    end

end