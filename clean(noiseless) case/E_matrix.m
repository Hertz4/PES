function [y,grad]=E_matrix(poles,GM_matrix,zM,N_orb)
NM=length(zM);
Num_poles=length(poles);
cvx_begin quiet
    variable R(N_orb,N_orb,Num_poles) hermitian semidefinite
    G3=[];
    for index=1:NM
        Ghere=zeros(N_orb);
        for m=1:Num_poles
            Ghere=Ghere+squeeze(R(:,:,m))/(1j*zM(index)-poles(m));
        end
        G3=[G3;Ghere];
    end
    obj=GM_matrix-G3;
    minimize( norm(obj,'fro'))
cvx_end
y=norm(obj,'fro');
if nargout > 1 % gradient required
    grad=zeros(1,Num_poles);
    for index=1:NM
        Ghere=obj(N_orb*(index-1)+1:N_orb*index,:);
        for k=1:Num_poles
            gradhere=squeeze(R(:,:,k))/((1j*zM(index)-poles(k))^2);
            grad(k)=grad(k)+real(sum(sum(conj(Ghere).*gradhere)));
        end
    end
    grad=-(1/y)*grad;
end

end