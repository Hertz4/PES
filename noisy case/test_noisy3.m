clear
% Code for PES method on noisy data, noise level = 1.6384e-2

%%
load("data3.mat")
% We have loaded the following quantities:
%load the Matsubara frequencies zM: a real-valued 1*NM array.
%   NM is the number of Matsubara frequencies.
%load the Matsubara data GM, which is NM*Norb*Norb complex-valued.
%   Norb is the number of orbitals
%   GM has noise level 1.6384e-2.
%load the real frequency mesh Omg, a real-valued 1*Nomg array.
%load the true value of Spectral function Spec_true, evaluated at Omg+1i*0.01.
%load the true value of Green's function Greenstrue, evaluated at Omg+1i*0.01.
%   Greenstrue is Nomg*Norb*Norb complex-valued.


fprintf("Noise level of Matsubara data is 1.6384e-2\n")

[NM,Norb,~] = size(GM);


% plot true spectral function
figure
plot(Omg,Spec_true,'linewidth',1.5)
title("True spectral function")
drawnow




%% Projection step 

% The projection here is done on the entire Norb*Norb matrix. 
% One can choose to do the projection only on the diagonal entries, which is faster.
% This grid could/should be changed to a finer one if needed.

x_grid = linspace(-6,6,200);
fprintf("Conducting Projection ...\n")
tic
[GM_proj, pol_proj, res_proj] = p_proj_matrix(GM, zM, x_grid,Norb);
t1 = toc;
fprintf(sprintf("Projection done in %f seconds!\n",t1))

%% Estimation step
eps_p=5e-4;
pol_ini=[];
for orb = 1:Norb
    [r,poles] = aaa(squeeze(GM_proj(:,orb,orb)),1j*zM);
    poles(abs(imag(poles))>eps_p)=[];
    pol = real(poles);
    pol_ini = [pol_ini;pol];
end
pol_ini = unique(sort(pol_ini));
size(pol_ini)
fprintf("Estimation of poles, done\n")



%% Semidefinite relaxation step
options = optimoptions('fminunc','Algorithm','quasi-Newton','SpecifyObjectiveGradient',true);
options.Display = 'iter';options.MaxIterations=300; 


fprintf("Conducting semidefinite relaxation fitting ...    ")
if options.Display == 'none'
    fprintf("The optimization details have been turned off.\n")
elseif options.Display == 'iter'
    fprintf("The optimization details are being displayed ...\n")
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

drawnow

subplot(2,1,2)
plot(Omg,real(squeeze(Greenstrue(:,1,3))),'linewidth',1.5)
hold on
plot(Omg,real(squeeze(Greens_calc(:,1,3))),'r--','linewidth',1.5)
legend('true value','calculation result')
title("Re(G_{13})")