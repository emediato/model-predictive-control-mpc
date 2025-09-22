clc
clear all

%% Parâmetros do conversor Buck
Vcc = 400;
Rd = 11.519;
Cd = 4.187e-6;
Cf = 4.187e-6;
L1 = 246.926e-6;
L2 = 246.926e-6;

%% Inicialização da função de transferência
s = tf('s');

%% Função de transferência do Filtro LCL + RC
Gvd = ((Vcc*(s*Rd*Cd+1))/(L1*L2*Cf*Cd*Rd))/(s*(s^3+(s^2)*(Cf+Cd)/(Rd*Cf*Cd) + s*(L1+L2)/(L1*L2*Cf) + (L1+L2)/(L1*L2*Cf*Cd*Rd)));

%% Diagrama de Bode do Filtro LCL + RC
figure;
bode(Gvd);
title('Diagrama de Bode - Função de Transferência de Tensão');
grid on;

%% ********************************
%% ********************************
%% ********************************
%% ********************************
%% ********************************

clc
clear all

%% Parâmetros do conversor Buck
Vp = 311; 
Vout = 380;  
R = 72.2;    
L = 1e-3; 
C = 50e-6; 
k_v = 1;
V_trip = 1;
k_pwm = 1;

%% Inicialização da função de transferência
s = tf('s');

%% Função de transferência de tensão
Gvd = (Vp/(2*Vout))*(R/(s*R*C +2));

%% Função de transferência de corrente 
Gid = Vout / (s*L);

%% Diagrama de Bode da planta de tensão
figure;
bode(k_pwm*k_v*Gvd);
title('Diagrama de Bode - Função de Transferência de Tensão');
grid on;

%% Diagrama de Bode da planta de corrente
figure;
bode(k_pwm*k_v*Gid);
title('Diagrama de Bode - Função de Transferência de Corrente');
grid on;

%% Especificações de projeto do controlador PI da tensão
Margem_Fase_Voltage = 120; 
Frequencia_Cruzamento_Voltage = 2.5; %%Hz
Frequencia_Cruzamento_Voltage_rad = 2*pi*Frequencia_Cruzamento_Voltage;

%% Cálculo do controlador PI da tensão
[mag, phase] = bode(k_pwm*k_v*Gvd, Frequencia_Cruzamento_Voltage_rad);

phase = phase(1); 
mag = mag(1); % Extrai o valor da magnitude

Angulo_PI_Voltage = Margem_Fase_Voltage - 180 - phase; 
w_z_Voltage = Frequencia_Cruzamento_Voltage_rad / tan(deg2rad(Margem_Fase_Voltage - 90 - phase));
k_c_Voltage = Frequencia_Cruzamento_Voltage_rad / (mag * sqrt(w_z_Voltage^2 + Frequencia_Cruzamento_Voltage_rad^2));

PI_Controller_Voltage = k_c_Voltage* (1 + s / w_z_Voltage);

display(PI_Controller_Voltage);

Constante_de_tempo_Voltage = 1/(w_z_Voltage);

%% Especificações de projeto do controlador PI da corrente
Margem_Fase_Current = 60; 
Frequencia_Cruzamento_Current = 2000; %%Hz
Frequencia_Cruzamento_Current_rad = 2*pi*Frequencia_Cruzamento_Current;

%% Cálculo do controlador PI da corrente
[mag, phase] = bode(k_pwm*k_v*Gid, Frequencia_Cruzamento_Current_rad);

phase = phase(1);
mag = mag(1); % Extrai o valor da magnitude

Angulo_PI_Current = Margem_Fase_Current - 180 - phase; 
w_z_Current = Frequencia_Cruzamento_Current_rad / tan(deg2rad(Margem_Fase_Current - 90 - phase));
k_c_Current = Frequencia_Cruzamento_Current_rad / (mag * sqrt(w_z_Current^2 + Frequencia_Cruzamento_Current_rad^2));

PI_Controller_Current = k_c_Current * (1 + s / w_z_Current);

display(PI_Controller_Current);

Constante_de_tempo_Current = 1/(w_z_Current);

%% ********************************
%% ********************************
%% ********************************
%% ********************************
%% ********************************
%% ********************************

%% ********************************

clear; clc;

%% PI CONTROLLER


 % Parameters
Vi = 400;           % Input voltage [V]
Ro = 31;            % Load resistance [Ohm]
Co = 10e-6;         % Output capacitance [F]
Lo = 1e-3;          % Output inductance [H]
kv = 1;             % Voltage gain
kpwm = 1;           % PWM gain
fs = 20e3;          % Switching frequency [Hz]

% Frequency range
w = logspace(1, 5, 1000); % 10 rad/s to 100 k rad/s

% Transfer function Gv(jw)
%s = 1j*w;
s = tf('s');

Gv = (Vi * 1/(Lo * Co)) / ( (s^2) + 1/(Ro * Co)*s + 1/(Lo * Co) );
bode(Gv);


%%
% Open-loop transfer function
FTLA_nc = kpwm .* Gv .* kv;

% Magnitude in dB
FTLA_nc_mod = 20 * log10(abs(FTLA_nc));

% Phase in degrees
FTLA_nc_fase = (180/pi) * angle(FTLA_nc);

% Plot Bode
figure;
subplot(2,1,1);
semilogx(w, FTLA_nc_mod);
grid on;
title('Magnitude of FTLA_{nc}');
ylabel('Magnitude (dB)');

subplot(2,1,2);
semilogx(w, FTLA_nc_fase);
grid on;
title('Phase of FTLA_{nc}');
ylabel('Phase (deg)');
xlabel('Frequency (rad/s)');




%% VSI OPEN LOOP
% Params
p = vsi_params();

% Frequency vector (rad/s): 1 Hz to 100 kHz
w = 2*pi*logspace(0, 5, 2000);



% Gv = Gv_vsi(w, p)
% Gv(jw): control-to-output transfer according to your formula
% w : rad/s (vector)
% p : struct from vsi_params()

s   = 1j * w(:);              % column vector
num = p.Vi * (p.Lo * p.Co);
den = (s.^2) + (p.Ro*p.Co).*s + (p.Lo*p.Co);
Gv  = num ./ den;             % column vector


%L = FTLA_nc_vsi(w, p)
% Open-loop transfer without controller: L = kpwm * kv * Gv(jw)
L  = p.kpwm * p.kv * Gv;

S = 1 ./ (1 + L);
T = L  ./ (1 + L);

% Plot closed-loop complementary sensitivity T (magnitude only)
figure;
semilogx(w, 20*log10(abs(T)), 'LineWidth', 1.4); grid on;
title('Closed-loop Complementary Sensitivity T(j\omega)');
xlabel('\omega (rad/s)');
ylabel('|T| (dB)');



% Plots magnitude (dB) and phase (deg) vs omega (rad/s) on log axis
% w  : rad/s vector
% H  : complex response (same size)
% ttl: title string (optional)

[mag_dB, phase_deg] = magphase(H);

figure;
subplot(2,1,1);
semilogx(w, mag_dB, 'LineWidth', 1.4); grid on;
ylabel('Magnitude (dB)');
title(ttl);

subplot(2,1,2);
semilogx(w, phase_deg, 'LineWidth', 1.4); grid on;
ylabel('Phase (deg)');
xlabel('\omega (rad/s)');


function p = vsi_params()
% Default parameters for VSI voltage loop
% From your spec:
%   Vi = 400 V, Ro = 31 ohm, Co = 10 uF, Lo = 1 mH
%   kv = 1, Vtrip = 1 V  -> kpwm = 1/Vtrip = 1
%   fs = 20 kHz
    p.Vi    = 400;          % [V]
    p.Ro    = 31;           % [ohm]
    p.Co    = 10e-6;        % [F]
    p.Lo    = 1e-3;         % [H]
    p.kv    = 1;            % [-]
    p.Vtrip = 1;            % [V]
    p.kpwm  = 1/p.Vtrip;    % [1/V]
    p.fs    = 20e3;         % [Hz]
end



