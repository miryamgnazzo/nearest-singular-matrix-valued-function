function [F_vec, F_der, Delta, iterations,h_vec, Norms,Diff_Norm] = inner_SP(Delta_in, A, f, mu, epsilon, tol, P, complex)
% Script for inner iteration, including the one for sparstity pattern
% It minimizes the functional, by using a constrained steepest descent method
% Delta_in = initial perturbations
% A = cell array containing the matrix coefficients
% f = function handling the analytic functions
% mu = points \mu we evaluate the functional at
% tol = tolerance we ask for
% P = cell array containing the structure of the perturbations
% complex = 0 if we want REAL distance
%         = 1 if we want COMPLEX distance

Norms = [];
Diff_Norm = [];

d = length(mu);
n = size(A{1},1); %size of the matrix

h  = epsilon/100; %h=stepsize
maxiter = 30000;
found = 0;
iterations = 0;
acceptedstep = 0;
F_vec = []; %vector for the evaluations F(eps)
h_vec = []; %vector for the stepsize h
h_vec = [h_vec; h];


Delta = Delta_in;
Delta = proj_SP(Delta, P);

Norm = norm(cell2mat(Delta), 'fro');

for i = 1: length(Delta)
    Delta{i} = Delta{i} / Norm;
end

while (iterations < maxiter) && (found == 0)

Sigma = zeros(1,d); %smallest singular values
U_fin = zeros(n,d); %associated left singular vectors
V_fin = zeros(n,d); %associated right singular vectors
    
    for k = 1:d
            
           M = zeros(n);
            
           eval_mu = f(mu(k));

           for t1 = 1: length(Delta)
              M = M + eval_mu(t1)*(A{t1}  + epsilon*Delta{t1});
           end
          
           [U,S,V] = svd(M);

           Sigma(k) = S(n,n);
           u = U(:,n);
           v = V(:,n);
           
           U_fin(:,k) = u;
           V_fin(:,k) = v;

    end
    
    F_epsilon = 0.5*sum((Sigma).^2);
    
    if (iterations == 0) 
       F_vec(1) = F_epsilon;
    end

        MM = repmat({zeros(n)},1, length(A));
    for j = 1: d
        eval_mu = f(mu(j));

        for r = 1: length(A)
            MM{r} = MM{r} + Sigma(j)*(eval_mu(r))'*U_fin(:,j)*V_fin(:,j)'; 
        end
    end

    MM_matrix = cell2mat(MM);
    
    if (complex == 0)
         MM_matrix = real(MM_matrix);
         MM = mat2cell(MM_matrix, n , n*ones(1,length(MM)));
    end

    MM = proj_SP(MM, P);

        MM_matrix = cell2mat(MM);

        F_der = - norm(MM_matrix, 'fro');
    
        Delta_matrix  = cell2mat(Delta); 

        eta = real(trace((Delta_matrix)'* MM_matrix));
        
        %STEP OF EULER METHOD
        Delta_matrix_new = Delta_matrix + h*(- MM_matrix + eta*Delta_matrix);

        Norm_pert = norm(Delta_matrix_new, 'fro');

        Delta_new = mat2cell(Delta_matrix_new, n, n*ones(1,length(A)));

        for y = 1: length(Delta)
            Delta_new{y} = Delta_new{y}/Norm_pert;
        end
               
        Sigma_new = zeros(1,d);
            
        for k = 1:d

           M_new = zeros(n);
            
           eval_mu = f(mu(k));

           for t2 = 1: length(Delta)
              M_new = M_new +  eval_mu(t2)*(A{t2}  + epsilon*Delta_new{t2});
           end

               [~,S_new,~] = svd(M_new);
               Sigma_new(k) = S_new(end,end);
        end
        
        F_epsilon_new = 0.5*sum((Sigma_new).^2);
    
   %CHECK on monotonicity of F_epsilon
    if ((F_epsilon_new) < (F_epsilon))
          
         %Accepted stepsize h
         h_vec = [h_vec; h];
         new_norm = norm(Delta_matrix_new - Delta_matrix ,'fro');
         Norms = [Norms; new_norm];
         Diff_Norm = [Diff_Norm; norm(Delta{1} - Delta_new{1},'fro')];
        
            if (new_norm <= tol) || (F_epsilon_new <= tol)
                found=1;
            end
    
            if (acceptedstep == 1) 
                 h = 1.2*h;
            end 
    
            acceptedstep = 1; 
            Delta = Delta_new;

    else
        acceptedstep = 0;
        h = (1/1.2)*h;
    end
   
   F_vec = [F_vec; F_epsilon_new];
   iterations = iterations + 1;
end
