function [poles,err,Rout]=SDR_FIT_boson_scalar(pol,G,zM,options)
func= @(pole) E_boson_scalar(pole,G,zM);
[poles,err] = fminunc(func, pol, options);
p_sign = sign(poles);
Num_poles=length(poles);
%NM = length(G);
Matrix = 1./(1j*zM - poles);

cvx_begin quiet
    variable R(1,Num_poles) nonnegative
    Rout = R.*p_sign';
    obj=Rout*Matrix-G;
    minimize( norm(obj,'fro'))
cvx_end

end