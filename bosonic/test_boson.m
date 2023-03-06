clear

load('cleandata.mat')
% We have loaded the following quantities:
%load the Matsubara frequencies zM: a real-valued 1*NM array.
%   NM is the number of Matsubara frequencies.
%load the clean Matsubara data GM, which is 1*NM complex-valued.
%load the real frequency mesh Omg, a real-valued 1*Nomg array.
%load the true value of Spectral function Spec_true, evaluated at Omg+1i*0.01.



%% add noise to GM
noise_level = 1e-6; % try changing to 0, 1e-6, 1e-4, 1e-2
GM = GM + max(abs(GM))*(randn(size(GM))+1j*randn(size(GM)))*noise_level;


%% Projection step
if noise_level~=0
    fprintf("Conducting projection ...\n")
    tic
    [GM_proj, poles, res] = p_proj_boson(transpose(GM), 1j*transpose(zM), ([-10:0.01:-0.01 0.01:0.01:10]'));
    t1 =toc;
    fprintf(sprintf("Projection step done in %f seconds!\n",t1))
else
    fprintf(sprintf("Projection step skipped because input data is clean!\n",t1))
end
%% Estimation step
eps_p=1;
if noise_level~=0
    [rat_approx, pol, res] = aaa(GM_proj, 1j*zM);
else
    [rat_approx, pol, res] = aaa(GM, 1j*zM);
end

pol(abs(imag(pol))>eps_p)=[]; pol(real(pol)==0)=[]; pol_ini = real(pol);
fprintf(sprintf("Estimation step done!\n",t1))

%% SDR fitting step


options = optimoptions('fminunc','Algorithm','quasi-Newton','SpecifyObjectiveGradient',true);
options.Display = 'none';options.MaxIterations=300;


fprintf("Conducting semidefinite relaxation fitting ...    ")
if options.Display == 'none'
    fprintf("The optimization details have been turned off.\n")
elseif options.Display == 'iter'
    fprintf("The optimization details are being displayed ...\n")
end
%pol_ini = pol_ini';

tic
[poles,err,Rout]=SDR_FIT_boson_scalar(pol_ini,GM,zM,options);
t1=toc;

fprintf(sprintf("SDR fitting done in %f seconds!\n",t1))

for i = 1:length(Omg)
    Spechere = 0;
    for m = 1:length(poles)
        Spechere = Spechere+Rout(m)/(Omg(i)+0.01*1j-poles(m));
    end
    Spec_calc(i) = -imag(Spechere)/pi;
end
    
figure
plot(Omg,Spec_real,'linewidth',2)
hold on
plot(Omg,Spec_calc,'r--','linewidth',3)
legend('True','Calc')
sgtitle(sprintf("noise-level = %.4e",noise_level))
drawnow