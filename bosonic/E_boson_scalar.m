function [y,grad]=E_boson_scalar(poles,G,zM)
    NM=length(zM);
    Num_poles=length(poles);
    p_sign = sign(poles);
    Matrix = 1./(1j*zM - poles);

    cvx_begin quiet
        variable R(1,Num_poles) nonnegative
        Rout = R.*p_sign';
        obj=G-Rout*Matrix;
        minimize( norm(obj,'fro'))
    cvx_end
    y=norm(obj,'fro');
    if nargout > 1 % gradient required
        grad=zeros(1,Num_poles);
        for index=1:NM
            Ghere=obj(index);
            gradhere=(R.*p_sign')./transpose((1j*zM(index)-poles).^2);
            grad=grad+real(conj(Ghere).*gradhere);
        end
        grad=-(1/y)*grad;
    end
  
end