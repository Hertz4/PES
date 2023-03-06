function [Fproj, poles, res] = p_proj_boson(F, Z, lam)
% Using convex optimization to project the data to the physical domain
% This is for the quasiparticle case


num_param = numel(lam);
A = zeros(numel(Z), numel(lam));
poles = lam;
for i = 1 : numel(Z)
  for j = 1 : num_param
    A(i,j) = 1.0/(Z(i)-poles(j));
  end
end

cvx_begin quiet
  variable res(num_param,1) nonnegative

  obj = A*(res.*sign(lam)) - F;
  sum(res.*sign(lam)) == 0;
  minimize(norm(obj, 'fro'))
cvx_end
Fproj = A*(res.*sign(lam));
end