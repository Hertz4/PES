clear

%%
load("data_clean.mat")
% We have loaded the following quantities:
%load the Matsubara frequencies zM: a real-valued 1*NM array.
%   NM is the number of Matsubara frequencies.
%load the Matsubara data GM, GM is NM*Norb*Norb complex-valued.
%   Norb is the number of orbitals
%load the real frequency mesh Omg, a real-valued 1*Nomg array.
%load the true value of Spectral function Spec_true, evaluated at Omg+1i*0.01.
%load the true value of Green's function Greenstrue, evaluated at Omg+1i*0.01.
%   Greenstrue is Nomg*Norb*Norb complex-valued.

[NM,Norb,~] = size(GM);


% plot true spectral function
figure
plot(Omg,Spec_true,'linewidth',1.5)
title("True spectral function")
drawnow




%% no Projection step since data is clean (noiseless)
fprintf("No Projection step since data is clean (noiseless).\n")
%% Estimation step
eps_p=1;
pol_ini=[];
GM_trace = zeros(size(squeeze(GM(:,1,1))));
for orb = 1:Norb
    GM_trace = GM_trace +squeeze(GM(:,orb,orb));
end
[r,poles] = aaa(GM_trace,1j*zM);
poles(abs(imag(poles))>eps_p)=[];
pol = real(poles);
pol_ini = [pol_ini;pol];
pol_ini = unique(sort(pol_ini));

fprintf("Estimation step: done.\n")





%% Semidefinite relaxation step
options = optimoptions('fminunc','Algorithm','quasi-Newton','SpecifyObjectiveGradient',true);
options.Display = 'iter';options.MaxIterations=300; 


fprintf("Conducting semidefinite relaxation fitting ...   ")
if options.Display == 'none'
    fprintf("The optimization details have been turned off. \n")
elseif options.Display == 'iter'
    fprintf("The optimization details are being displayed.\n")
end
tic
[polM_here,XM,errM, Spec_calc,Greens_calc] = SDR_FIT(pol_ini,GM,zM,Norb,Omg,options);
t1 = toc;
fprintf(sprintf("SDR fitting done in %f seconds!\n",t1))
    

%% plot the result for spectral functions and non-diagonal entry G13
figure
subplot(2,1,1)
plot(Omg,Spec_true,'linewidth',1.5)
hold on
plot(Omg,Spec_calc,'r--','linewidth',1.5)
legend('true spectrum','calculation result')
title("Spectrum")


subplot(2,1,2)
plot(Omg,real(squeeze(Greenstrue(:,1,3))),'linewidth',1.5)
hold on
plot(Omg,real(squeeze(Greens_calc(:,1,3))),'r--','linewidth',1.5)

legend('true value','calculation result')
title("Re(G_{13})")