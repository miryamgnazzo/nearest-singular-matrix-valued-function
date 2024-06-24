function P = sparsity_pattern(A)
%% Recognize sparsity pattern and collect it in the 1-0 cell array P
% A = cell array of the coefficient

N = size(A{1},1);

AA = cell2mat(A);
[n1,n2] = size(AA);

PP = zeros(n1,n2);

for i = 1 : n1
   for  j = 1 : n2
      if (AA(i,j) ~= 0)
         PP(i,j) = 1;
      end
   end
end

P = mat2cell(PP, N, N*ones(1, length(A)) );