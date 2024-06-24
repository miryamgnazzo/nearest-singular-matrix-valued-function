function M_Proj = proj_SP(M, P)
%% Project onto the structure (sparsity pattern)
% P = Cell array containing sparsity structure (1-0 structure)
% M = Cell array containing the coefficient matrices to project

n = size(M{1},1);

MM = cell2mat(M);
PP = cell2mat(P);

[n1,n2] = size(MM);

MM_Proj = MM;

for i=1:n1
   for j=1:n2
      if (PP(i,j)==0)
          MM_Proj(i,j)=0;
      end
   end
end

M_Proj = mat2cell(MM_Proj, n, n*ones(1,length(M)));

