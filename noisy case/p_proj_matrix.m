function [Fproj, poles, res] = p_proj_matrix(F, Z, lam,dim)
% Using convex optimization to project the data to the physical domain
% This is for the particle case, matrix case
%F: Nw*dim*dim
Z = 1j*reshape(Z,[],1);
num_param = numel(lam);
A = zeros(numel(Z), numel(lam));
poles = lam;
for i = 1 : numel(Z)
	for j = 1 : num_param
        A(i,j) = 1.0/(Z(i)-poles(j));
	end
end



cvx_begin quiet
	variable res(dim,dim,num_param) hermitian semidefinite
    for j1 = 1:dim
        for j2 = 1:dim
            Fproj(:,j1,j2) = A*squeeze(res(j1,j2,:)) ;
        end
    end
    sum(res,3) == eye(dim);
    obj = Fproj-F;
    minimize(norm(obj(:), 'fro'))
cvx_end
for j1 = 1:dim
    for j2 = 1:dim
        Fproj(:,j1,j2) = A*squeeze(res(j1,j2,:)) ;
    end
end
end