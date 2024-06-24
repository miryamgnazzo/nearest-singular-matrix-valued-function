function [A_new, DD, Max_det] = normalization_det(A, f, N)
%% Normalization of the determinant on the unit disk
% A = cellarray with the coefficients
% f = analytic functions
% N = number of points

mu = exp((2*pi*1i).*((1:N)./N));

s = size(A{1},1);

vec_det = zeros(1,N);

for i = 1 : N
    
   eval_mu = f(mu(i));

   Matrix = zeros(s,s);

   for j = 1 : length(eval_mu)
        Matrix = Matrix + eval_mu(j)*A{j};
   end
   
   [U,S,V] = svd(Matrix);
   vec_det(i) = det(U)*det(V)*(prod(diag(S)));
end

Max_det = abs(max(vec_det));

DD = nthroot(Max_det,s);

A_new = A;

for k = 1 : length(A_new)
    A_new{k}  = A_new{k}*(1/DD);
end