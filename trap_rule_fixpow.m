function [approx, val] = trap_rule_fixpow(A, f, exponent)
%% Application of the trapezoidal rule for the approximation of the integral
% A = cell array containing the matrix coefficients
% f = function handling the analytic functions
% exponent = exponent we are considering

maxiter = 7;
N = 20;
val = [];

z = exp((2*pi*1i/N)*(1:N));

eval = evaluation_det(A, f, z);

powers = z.^(-exponent);
approx = sum(powers.*eval);
approx = abs(approx/N);

k = 2;
while (k <= maxiter)
    N = N*k;
    z = exp((2*pi*1i/N)*(1:N));
    eval = evaluation_det(A, f, z);
    
    powers = z.^(-exponent);
    approx_new = sum(powers.*eval);
    approx_new = abs(approx_new/N);
    
    if (abs(approx-approx_new) <= 10^-11)
        k = 20;
    end

    approx = approx_new;
    val = [val approx];
    k = k+1;
end
