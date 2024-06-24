%EXAMPLE 

n = 2;
rng(2)

A2 = randn(n);
A1 = randn(n);
A0 = randn(n);

A = {A2, A1, A0};

threshold = 1e-3;
tol_appr = 1e-12;
tol_sing=1e-10;

complex = 0;
str = 0;

[distance,Pert,counter,iter_tot,fin_points,max_needed_points,DD]=struct_distance(A, @f, threshold, tol_appr, tol_sing, str, complex);

function fv = f(x)
 fv = [x, exp(-x), 1];
end
