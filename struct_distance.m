function [distance,Pert,counter,iter_tot,fin_points,max_needed_points,DD]=struct_distance(A, f, threshold, tol_appr, tol_sing, str, complex)
% APPROXIMATION DISTANCE TO SINGULARITY FOR MATRIX-VALUED FUNCTIONS
% A = cell array containing the matrix coefficients
% f = function handling the analytic functions
% threshold = threshold for recompute the number of points
% tol_appr = tolerance we ask in the approximation of the determinant
% tol_sing = tolerance we ask on each singular value
% str = 0 if we want the UNSTRUCTURED distance to singularity
%     = 1 if we want the STRUCTURED distance to singularity (SPARSITY)
% complex = 0 if we want REAL distance
%         = 1 if we want COMPLEX distance

n = size(A{1},1);

%NORMALIZATION of the initial coefficients
[A, DD, ~] = normalization_det(A, f, 20);

P = repmat({zeros(n)},1,length(A));
 
%UNSTRUCT vs STRUCT distance
    if (str == 0)
        for i1 = 1 : length(A)
            P{i1} = ones(n);
        end
    else
        P = sparsity_pattern(A);
    end


%FIND the initial number of needed points
[d , coef_Tay] = find_exp(A, f, tol_appr);

max_needed_points = d;

fprintf('number of STARTING points %d with |a_n| %d \n \n ', d, coef_Tay)

mu = exp((2*pi*1i).*((1:d)./d));

%Initialization of [Delta2, Delta1, Delta0]
Sigma = zeros(1,d); %smallest singular values
U_fin = zeros(n,d); %associated left singular vectors
V_fin = zeros(n,d); %associated right singular vectors

    for k = 1:d

          eval_mu = f(mu(k));

           M = zeros(n);
           for y = 1 : length(A)
                M = M + eval_mu(y)*A{y};
           end
    
           [U,S,V] = svd(M);
           Sigma(k) = S(end,end);

           u = U(:,end);
           v = V(:,end);

           U_fin(:,k) = u;
           V_fin(:,k) = v;

    end
    

MM = repmat({zeros(n)},1,length(A));

for j = 1 : d
    
    eval_mu = f(mu(j));

    for r = 1 : length(A)
        MM{r} = MM{r} + Sigma(j)*(eval_mu(r))'*U_fin(:,j)*V_fin(:,j)';
    end

end  

MM_matrix = cell2mat(MM);

if (complex == 0)
     MM_matrix = real(MM_matrix);
     MM = mat2cell(MM_matrix, n , n*ones(1,length(MM)));
end

if (str == 1)
        MM = proj_SP(MM ,P);
end

MM_matrix = cell2mat(MM);

NN = norm(MM_matrix,'fro');

Delta = MM;

for t = 1: length(MM)
    Delta{t} = (-1/NN)*Delta{t};
end

epsilon = 10^-3;

%Choice of epsilon
F_epsilon = 0.5*sum(((Sigma)).^2);
F_der_epsilon = - NN;

epsilon = -(F_epsilon)/(F_der_epsilon);


Norm_A = norm(cell2mat(A), 'fro');

if (epsilon >= Norm_A)
    epsilon = 10^-5;
end

value_epsilon = [];
value_epsilon = [value_epsilon epsilon];

%choice of tolerances
tol_inner = max((F_epsilon)/1e3,1e-9);
tol_out= tol_sing*d;

%INNER ITERATION - FIRST ONE
fprintf('INNER Tollerance %d \n', tol_inner)

[F_vec, F_der, Delta, ~,~, ~,~ ] = inner_SP(Delta, A, f, mu, epsilon, tol_inner, P, complex);

Eval_F_eps = []; %Vector for the G(eps_k)
Distanze = [];  %Vector for the distance between ep_up- ep_low at each step k

%Inizializzo estremi dell'intervallo [ep_low, ep_up]
ep_low = epsilon;
ep_up = Norm_A;

eps_final = epsilon;     
iter_tot = 0;
dist_epsilon = 1;
counter=0;

while (dist_epsilon > 1e-6)&&(iter_tot<=30)
    
    tol_inner=max((F_epsilon)/(10^3),1e-10);
   
    M_old = [];

    %STEP k of the method
    for c1 = 1: length(A) 
        M_old{c1} = A{c1} + epsilon*Delta{c1};
    end
    
    %INNER ITERATION 
    fprintf('INNER Tollerance %d \n', tol_inner)

    [F_vec, F_der, Delta, ~,~, ~,~] = inner_SP(Delta, A, f, mu, epsilon, tol_inner, P, complex);

    %STEP k+1 of the method
    M_new = [];

    %STEP k of the method
    for c2 = 1: length(A) 
        M_new{c2} = A{c2} + epsilon*Delta{c2};
    end
        
    F_epsilon = F_vec(end);
    
    Eval_F_eps = [Eval_F_eps, F_epsilon];

    if (F_epsilon > tol_out)
        ep_low = epsilon;
        epsilon_new = epsilon - ((F_epsilon - tol_out)/(F_der));
    else
       ep_up = epsilon;
       epsilon_new = (ep_low+ep_up)/2;
    end

    dist_epsilon = abs(epsilon_new -epsilon);
    
   %New epsilon - Computation    
    if (epsilon_new > ep_low) && (epsilon_new < ep_up)
        epsilon = epsilon_new;
    else
       epsilon = (ep_low+ep_up)/2;
    end
    
    
    %Check on the increments
    Diff = cell2mat(M_old) - cell2mat(M_new);
    
    increase= norm(Diff, 'fro');
    increase = increase/norm(cell2mat(M_old),'fro');
    
    if (increase > threshold)
        
        %CHECK on the number of needed points
        [d_check , coef_check] = find_exp(M_new, f, tol_appr);
        
        fprintf('INCREASE %d ---> \n', increase)
        fprintf('number of NEW points %d with |a_n| %d \n \n ', d_check, coef_check)
        
        if (d_check >= 2*d)
            counter = counter+1;
            d = 2*d;            
            fprintf('now we have NEW points %d \n \n ', d)
            
            mu = exp((2*pi*1i).*((1:d)./d));
            tol_out = tol_sing*d
        else
            if (d_check <= 0.5*d)
                counter = counter+1;
                d = floor(0.5*d);
                fprintf('now we have NEW points %d \n \n ', d)
                mu = exp((2*pi*1i).*((1:d)./d));
                tol_out = tol_sing*d
            end
        end
        
        %needed points
        if (d > max_needed_points)
            max_needed_points = d;
        end
        
    end
    
    value_epsilon = [value_epsilon epsilon];
    Distanze = [Distanze, dist_epsilon];

   iter_tot = iter_tot+1;
end

fprintf('------------- \n')
fprintf('NUMBER OF RE_COMPUTATIONS %d on ITER %d:  \n \n', counter, iter_tot)

eps_final = epsilon;

%FINAL RESULTS
distance = eps_final;
Pert = Delta;
fin_points = d;
