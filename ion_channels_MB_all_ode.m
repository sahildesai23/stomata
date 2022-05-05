function dydt = ion_channels_MB_all_ode(t,y, p, f_L)

Vm    = y(1);
Clin  = y(2); 
Kin   = y(3);
L     = f_L(t);

%%
Cm  = p.Cm;  % Membrane capacitance
F   = p.F;   % Faraday constant
vol_in = p.vol_in; % Guard cell volume
area = p.area;
R = p.R;
T = p.T;

Km = p.Km;
h = p.h;
P_Cl_max = p.P_Cl_max; 
P_K_max = p.P_K_max; % m/s


NSLAC1 = p.NSLAC1;
NGORK  = p.NGORK;

V1_2 = p.V1_2;
S    = p.S;
%%
p_open_Cl = (L^h)/(Km^h + L^h);

p_open_K  = 1/(1 + exp((V1_2 - Vm)/S)).^2;


%%

I_Cl_prefactor = NSLAC1*p_open_Cl*(P_Cl_max*Vm*F/(R*T));
I_Cl_num = 1e3*(Clin - p.Clout_0*exp(Vm*F/(R*T)));
I_Cl_den = 1 - exp(Vm*F/(R*T));

I_Cl_GHK = I_Cl_prefactor*I_Cl_num/I_Cl_den; % 1e3 multiplied to convert mol/L to mol/m^3

% I_Cl_GHK is in A/m^2

%%
I_K_prefactor = NGORK*p_open_K*(P_K_max*Vm*F/(R*T));
I_K_num       = 1e3*(Kin - p.Kout_0*exp(-Vm*F/(R*T)));
I_K_den       = 1 - exp(-Vm*F/(R*T));

I_K_GHK = I_K_prefactor*I_K_num/I_K_den; % 1e3 multiplied to convert mol/L to mol/m^3

% I_K_GHK is in A/m^2

%%
dVmdt   = -1*(1/Cm)*(I_Cl_GHK + I_K_GHK)*area;
dClindt = (I_Cl_GHK*area)/(F*vol_in);
dKindt  = -(I_K_GHK*area)/(F*vol_in);

dydt = [dVmdt; dClindt; dKindt];

end