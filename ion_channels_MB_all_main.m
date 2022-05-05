clc; clear; close all;
%%
% Reaction parameters
p.Km = 5.2; % muM
p.h = 1; % Hill Coefficient
% p.GCl_max = 135e-12;  % ohm^-1
p.NSLAC1 = 100;
p.NGORK  = 100;


% membrane capacitance and reversal potentials
p.Cm = 6e-12; % Farads, Membrane capacitance Schroeder 1987
p.F = 96485.3329; % Coulombs/mol
p.vol_in = 1*1e-18*1e3; % L
p.vol_out = 1000*1e-18*1e3;
p.R = 8.317;
p.T = 298;
p.area = 1009e-12;



p.Clin_0  = 80e-3; % M
p.Clout_0 = 25e-3; % M

p.Kin_0  = 150e-3; % M
p.Kout_0 = 50e-3; % M

p.P_Cl_max = 2.45e-6; % m/s 
p.P_K_max  = 15e-6; % m/s

% GORK parameters
p.S = 21e-3;
p.V1_2 = -7e-3;


% Time scales
tau_electric  = p.R*p.T*p.Cm/(p.Clout_0*p.area*p.P_Cl_max*p.F^2)
tau_diffusion = 1e-3*p.vol_in/(p.area*p.P_Cl_max)

%%
% OST1* ligand concentration changes

Li = 0;
Lf = 500; % muM
tstep = 1; 
steepness = 0;
f_L = @(t)(Li + Lf*(t>tstep).*((t-tstep).^2)./(steepness+(t-tstep).^2));

Vm0 = -150e-3; % V
IC = [Vm0 p.Clin_0 p.Kin_0];  % [Vm, n]
tspan = [0 3]; 

[t,y] = ode15s(@(t,y) ion_channels_MB_all_ode(t,y,p, f_L),tspan, IC);

%%
Vm    = y(:,1);
Clin  = y(:,2);
Kin    = y(:,3);

%%

L = f_L(t);   

p_open_Cl = (L.^p.h)./(p.Km^p.h + L.^p.h);

%gCl = p.GCl_max*p_open_Cl;
ECl = -0.025*log(p.Clout_0./Clin); % V
EK  =  0.025*log(p.Kout_0./Kin); % V

I_Cl_prefactor = p.NSLAC1*p_open_Cl.*(p.P_Cl_max*Vm*p.F/(p.R*p.T));
I_Cl_num = Clin - p.Clout_0*exp(Vm*p.F/(p.R*p.T));
I_Cl_den = 1 - exp(Vm*p.F/(p.R*p.T));

ICl = 1e3*(I_Cl_prefactor.*I_Cl_num./I_Cl_den)*p.area;

%%
p_open_K  = 1./(1 + exp((p.V1_2 - Vm)/p.S)).^2;


I_K_prefactor = p.NGORK*p_open_K.*(p.P_K_max*Vm*p.F/(p.R*p.T));
I_K_num = Kin - p.Kout_0*exp(-Vm*p.F/(p.R*p.T));
I_K_den = 1 - exp(-Vm*p.F/(p.R*p.T));

IK = 1e3*(I_K_prefactor.*I_K_num./I_K_den)*p.area;




Itot = ICl + IK;



deltaCtot = (Clin) - (p.Clin_0) + (Kin) - p.Kin_0;

% 1mM = 0.0024776 MPa
deltapi = deltaCtot*1e3*0.0024776; % MPa

%
dClindt = ICl/(p.F*p.vol_in);



%%
figure
subplot(10,1,1)
plot(t,p_open_Cl,'-', 'Linewidth', 3)
xlabel('time (s)')
ylabel('P_{0}^{Cl}')
set(gca,'FontSize',14)

subplot(10,1,2)
plot(t,Vm*1e3,'-', 'Linewidth', 3)
xlabel('time (s)')
ylabel('V_m (mV)')
set(gca,'FontSize',14)

subplot(10,1,3)
plot(t,Clin*1e3,'-', 'Linewidth', 3)
xlabel('time (s)')
ylabel('[Cl^-]_{in} (mM)')
set(gca,'FontSize',14)

subplot(10,1,4)
plot(t,ICl*1e12,'-', 'Linewidth', 3)
xlabel('time (s)')
ylabel('I_{Cl} (pA)')
set(gca,'FontSize',14)

subplot(10,1,5)
plot(t,p_open_K,'-', 'Linewidth', 3)
xlabel('time (s)')
ylabel('P_{0}^{K}')
set(gca,'FontSize',14)

subplot(10,1,6)
plot(t,Kin*1e3,'-', 'Linewidth', 3)
xlabel('time (s)')
ylabel('[K^+]_{in} (mM)')
set(gca,'FontSize',14)

subplot(10,1,7)
plot(t,IK*1e12,'-', 'Linewidth', 3)
xlabel('time (s)')
ylabel('I_{K} (pA)')
set(gca,'FontSize',14)

subplot(10,1,8)
plot(t,Itot*1e12,'-', 'Linewidth', 3)
xlabel('time (s)')
ylabel('I_{Tot} (pA)')
set(gca,'FontSize',14)


subplot(10,1,9)
plot(t,ECl*1e3,'-', 'Linewidth', 3)
xlabel('time (s)')
ylabel('E_{Cl} (mV)')
set(gca,'FontSize',14)

subplot(10,1,10)
plot(t,EK*1e3,'-', 'Linewidth', 3)
xlabel('time (s)')
ylabel('E_{K} (mV)')
set(gca,'FontSize',14)
%%


figure
plot(Vm*1e3,Clin*1e3,'-', 'Linewidth', 3)
xlabel('Vm (mV)')
ylabel('[Cl^-]_{in} (mM)')
set(gca,'FontSize',14)

%%

figure
plot(t,deltapi,'-', 'Linewidth', 3)
xlabel('time (s)')
ylabel('\Delta\Pi (MPa)')
set(gca,'FontSize',14)


%%

figure
plot(t,dClindt*1e3 ,'-', 'Linewidth', 3)
xlabel('time (s)')
ylabel('d[Cl^-]/dt (mM/s)')
set(gca,'FontSize',14)

%%

figure
plot(Vm*1e3,dClindt*1e3 ,'-', 'Linewidth', 3)
xlabel('Vm (mV)')
ylabel('d[Cl^-]/dt (mM/s)')
set(gca,'FontSize',14)


