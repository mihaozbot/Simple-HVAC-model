
load('data_HVAC')

HVAC.ts = 1; %Sampling time

HVAC.k_z = 1; 20; %Zone temperature parameter
HVAC.c_z = 1; 12e4;

%Disturbance temperature
HVAC.K_D = 0.1;
HVAC.tau_D = 100;
HVAC.d_D = 0;

%Input temperature
HVAC.K_I = 0.02;
HVAC.tau_I = 100;
HVAC.d_I = 3;
HVAC.a_I = 1;
HVAC.b_I = 2;

%Outside temperature
HVAC.T_o = @(k) 295 + 10*sin(k*(2*pi)/1000) + 30*sin(k*(2*pi)/10000);

HVAC.Q_I = @(Q_I, T_z, T_O_d_I, u) Q_I + HVAC.ts*(HVAC.K_I*((HVAC.a_I+HVAC.b_I*(T_O_d_I-T_z)))*HVAC.k_z*atan(u)-Q_I)/HVAC.tau_I; %Energy input
HVAC.Q_D = @(Q_D, T_z, T_O_d_D) Q_D + HVAC.ts*(HVAC.K_D*(T_O_d_D-T_z)-Q_D)/HVAC.tau_D; %Energy disturbance
HVAC.T_z = @(T_z, Q_I, Q_D) T_z + HVAC.ts*(Q_I/HVAC.c_z + Q_D/HVAC.c_z); %Zone temperature

