function eval = evaluation_det(A,f, vec_lam)
% Evaluation of the determinant at points vec_lam
% A = cell array containing the matrix coefficients
% f = function handling the analytic functions
% vec_lam = points we evaluate at

s = length(vec_lam);
eval = zeros(1,s);

t = size(A{1},1);

for i = 1:s

    lam = vec_lam(i);
    eval_lam = f(lam);
    
     Matrix = zeros(t);

        for j = 1 : length(A)
            Matrix = Matrix + eval_lam(j)*A{j};
        end
    
    [U,S,V] = svd(Matrix);
    eval(i) = det(U)*det(V')*prod(diag(S));
end
