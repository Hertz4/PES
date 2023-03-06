function [polM_here,XM,errM, Spec_calc,Greens_calc] = SDR_FIT(pol_ini,GM,zM,Norb,Omg,options)

Gw4d = [];
for k = 1:length(zM)
    Gw4d = [Gw4d;squeeze(GM(k,:,:))];
end

%[polM_here,XM,errM] = Hyb_fit_matrix(pol_ini,Gw4d,zM,Norb,options);

func= @(pole) E_matrix(pole,Gw4d,zM,Norb);
[polM_here,errM] = fminunc(func, pol_ini, options);


NM=length(zM);
Num_pol=length(polM_here);
cvx_begin quiet
    variable XM(Norb,Norb,Num_pol) hermitian semidefinite
    G3=[];
    for index=1:NM
        Ghere=zeros(Norb);
        for m=1:Num_pol
            Ghere=Ghere+squeeze(XM(:,:,m))/(1j*zM(index)-polM_here(m));
        end
        G3=[G3;Ghere];
    end
    obj=Gw4d-G3;
    minimize(norm(obj,'fro'))
cvx_end
y=norm(obj,'fro');

Spec_calc = zeros(1,length(Omg));
Greens_calc = zeros(length(Omg),Norb,Norb);
for i = 1 : length(Omg)
    Greenhere = zeros(Norb,Norb);
    for l = 1: length(polM_here)
        Greenhere = Greenhere + (1.0/(Omg(i)+0.01*1j-polM_here(l)))*squeeze(XM(:,:,l));
    end
    Spec_calc(i) = -imag(trace(Greenhere))/pi;
    Greens_calc(i,:,:) = Greenhere;
end