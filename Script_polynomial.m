%EXAMPLE - POLYNOMIAL
%It requires the nlevp package

[A, f] = nlevp('damped_beam', 20);

threshold = 1e-3;
tol_appr = 1e-12;
tol_sing=1e-10;

complex = 0;
str = 1;

[distance,Pert,counter,iter_tot,fin_points,max_needed_points,DD]=struct_distance(A, f, threshold, tol_appr, tol_sing, str, complex);
