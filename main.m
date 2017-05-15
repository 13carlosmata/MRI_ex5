%% Simulation Exercise number5
clear all
close all
addpath 'seemri'
%% Changing TEs
IV = waterandfat(3,2);
tau = 2e-3;
TE_f=[5e-3, 50e-3, 100e-3];
S1={};
k1={};
figure
for i=1:length(TE_f)
    [S,kmax,~] = ex4(TE_f(i),5,tau,IV);
    S1{i}=S;
    k1{i}=kmax;
    fprintf("Done- \n")
    subplot(1,3,i)
    mrireconstruct(S1{i},k1{i},'Plot',true);
    title(['TE = ',num2str(TE_f(i)), 's'])
end
suptitle('Contrast changes using brain2mmpixelfuzzy at different TEs');
%% Changing TRs
TR_f=[2,1,0.5];
S2={};
k2={};
figure
for i=1:length(TR_f)
    [S,kmax,~] = ex4(5e-3,TR_f(i),tau,IV);
    S2{i}=S;
    k2{i}=kmax;
    fprintf("Done- \n")
    subplot(1,3,i)
    mrireconstruct(S2{i},k2{i},'Plot',true);
    title(['TR = ',num2str(TR_f(i)), 's'])
end
suptitle('Contrast changes using brain2mmpixelfuzzy at different TRs');
%% Computing TR and TE
T1_gm = 833e-3; T2_gm = 83e-3; p_gm = 0.86;
T1_wm = 500e-3; T2_wm = 70e-3; p_wm = 0.77;
D1=[]; D2=[];
T=0.001:0.001:6;
for i=1:length(T)
    I_gm_1 = p_gm*(1-exp(-T(i)/T1_gm));
    I_wm_1 = p_wm*(1-exp(-T(i)/T1_wm));
    d_1 = abs(I_gm_1 - I_wm_1);
    I_gm_2 = p_gm*exp(-T(i)/T2_gm);
    I_wm_2 = p_wm*exp(-T(i)/T2_wm);
    d_2 = abs(I_gm_2 - I_wm_2);
    D2 = [D2;d_2];
    D1 = [D1;d_1];
end
[y,i]=max(D1);
fprintf (['For a T1-weighted image, the optimal TR value is at: ',num2str(T(i)),'s \n'])
[y,j]=max(D2);
fprintf (['For a T2-weighted image, the optimal TE value is at: ',num2str(T(j)),'s \n'])
%% Graphing for a T1-weighted image and T2-weighted image
% - T1-weighted image
TR_t1 = T(i);
TE_t1 = TR_t1/10;
figure
subplot(1,2,1)
[S_t1,kmax_t1,~] = ex4(TE_t1,TR_t1,tau,IV);
mrireconstruct(S_t1,kmax_t1,'Plot',true);
title(['T1-weighted image at TE:',num2str(TE_t1),'s and TR:',num2str(TR_t1),'s'])
% - T2-weighted image
TE_t2 = T(j);
TR_t2 = TE_t2*10;
subplot(1,2,2)
[S_t2,kmax_t2,~] = ex4(TE_t2,TR_t2,tau,IV);
mrireconstruct(S_t2,kmax_t2,'Plot',true);
title(['T2-weighted image at TE:',num2str(TE_t2),'s and TR:',num2str(TR_t2),'s'])
suptitle('Q3. Using Optimal Parameters')
%% Chemical Shift
tau_shift = [0.002,0.004,0.01];
tp = 1e-3;
TE_shift = 1.5*(2*tau_shift+tp);
TR_shift = 0.5;
figure
for i=1:length(tau_shift)
    [S_s,kmax_s,S_iv] = ex4(TE_shift(i),TR_shift,tau_shift(i),IV);
    subplot(1,3,i)
    mrireconstruct(S_iv,kmax_s,'Plot',true);
    title(['tau = ',num2str(tau_shift(i)), 's'])
end
suptitle('Chemical shift artifact using the waterandfat phantom at different tau values')
%%  Q5 Optional
gammabar = 42.58e6;  gamma = 2*pi*gammabar;
dir = 3.5e-6; B0 = 1.5;
delta_wc = gamma * dir * B0;
delta_x = 0.001;
G_new = dir*B0/delta_x;
fprintf(['For a chemical shift lower than 1mm, a G greater than ',num2str(G_new),'T/mm must be used \n']);
%% Assignment 4 code
function [S,kmax,S_iv] = ex4(TE,TR,tau,iv)
addpath 'seemri'
%%  Preparations
gammabar = 42.58e6;  gamma = 2*pi*gammabar;
%% Basic Gradient Echo Imaging
res = 2;
B0 = 1.5; tp =  1e-3; alpha = pi/2; B1 = alpha / (gamma*tp); f_rf = gammabar*B0;
rf = RectPulse(B1, f_rf, 0, tp);
% iv = disc(3,res);
% TE = 5e-3;
% TR = 5;
% tau = 2e-3;
FOV= 220;   %Changed to 220 when using brain* functions    
dw = res;
kmax =  1/(2*dw);
dk = 1/FOV;
dk = kmax/ceil(kmax/dk);
ks = -kmax:dk:kmax-dk;
w = gamma*B0;
Gpexs = ks/(gammabar*tau);
gx = Gradient([tp tp+tau], {Gpexs 0});
Gfey1 = -kmax/(gammabar*tau);
Gfey2 = kmax/(gammabar*tau);
gy = Gradient([tp tp+tau TE-tau TE+tau], [Gfey1 0 Gfey2 0]);
%%
dt = dk/(gammabar*Gfey2);
adc = ADC(TE-tau, TE+tau, dt);
[S_iv,t_iv] = seemri(iv,B0,rf,gx,gy,adc,TR,length(Gpexs),'Plot',false);
[S,ts] = brain_2mm_pixel_fuzzy(B0,rf,gx,gy,adc,TR,length(Gpexs));
end
