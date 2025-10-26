% Parameters from the first image
V_in = 48;          % Bus voltage [V]
R_a = 0.016;        % Armature resistance [Ω]
L_a = 0.06;         % Armature inductance [H]
V_trip = 1;         % Triangular peak value [V]

% Controller gains
k_i = 1;            % Current sensor gain
kia = 1;            % Current controller gain (?)
kpwm = 1/V_trip;    % PWM gain

% Transfer function G_id(s) = i_a(s)/d(s) = V_in / (R_a + s*L_a)
% Create the transfer function
numerator = V_in;
denominator = [L_a, R_a];
G_id = tf(numerator, denominator);

% Open-loop transfer function FTLA
FTLA_net = k_i * kpwm * G_id;

% Create frequency vector (adjust range as needed)
omega = logspace(0, 5, 1000); % 10^0 to 10^5 rad/s

% Calculate magnitude and phase
[mag, phase] = bode(FTLA_net, omega);
mag = squeeze(mag);
phase = squeeze(phase);

% Convert to dB for magnitude plot
mag_dB = 20*log10(mag);

% Create Bode plot
figure;

% Magnitude plot
subplot(2,1,1);
semilogx(omega, mag_dB, 'b', 'LineWidth', 2);
grid on;
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
title('Bode Diagram - Magnitude');
legend('FTLA_{net}');

% Phase plot
subplot(2,1,2);
semilogx(omega, phase, 'r', 'LineWidth', 2);
grid on;
xlabel('Frequency (rad/s)');
ylabel('Phase (degrees)');
title('Bode Diagram - Phase');
legend('FTLA_{net}');

% Alternative: Use MATLAB's built-in bode plot
figure;
bode(FTLA_net);
grid on;
title('Bode Diagram using MATLAB bode function');

% Display transfer function information
disp('Transfer Function G_id(s):');
disp(G_id);
disp('Open-loop Transfer Function FTLA_net(s):');
disp(FTLA_net);
disp(['DC Gain: ', num2str(dcgain(FTLA_net))]);
disp(['Crossover Frequency: ', num2str(bandwidth(FTLA_net)), ' rad/s']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
% closed loop current
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parâmetros do Sistema  
V_in = 48;              % Tensão de barramento [V]
R_a = 0.016;            % Resistência de armadura [Ω]
L_a = 0.06;             % Indutância de armadura [H]
J = 0.018;              % Momento de inércia [kg·m²]
L_af = 0.078;           % Indutância mútua [H]
I_f = 1.6;              % Corrente de campo [A]
V_trip = 1;             % Valor de pico da triangular [V]

% Ganhos do controlador
k_i = 1;                % Ganho do sensor de corrente
kia = 1;                % Ganho do controlador de corrente
kpwm = 1/V_trip;        % Ganho do PWM

%% 1. MALHA DE CORRENTE - Projeto do Controlador PI
fprintf('=== PROJETO DO CONTROLADOR PI PARA MALHA DE CORRENTE ===\n');

% Frequências de projeto
f_i = 1.66e4;           % Frequência do conversor [Hz]
f_u = f_i / 10;         % Frequência de cruzamento [Hz]
w_d = 2 * pi * f_u;     % Frequência angular de cruzamento [rad/s]

fprintf('Frequência de cruzamento: f_u = %.2f Hz\n', f_u);
fprintf('Frequência angular de cruzamento: w_d = %.2f rad/s\n', w_d);

% Margem de fase desejada
M_p = 60 * pi/180;      % 60 graus em radianos
fprintf('Margem de fase desejada: M_p = %.2f graus\n', M_p*180/pi);

% Função de transferência da planta de corrente
G_id = tf(V_in, [L_a, R_a]);
fprintf('\nFunção de transferência G_id(s):\n');
disp(G_id);

% FTLA inicial (malha aberta sem compensador)
FTLA_int = k_i * kpwm * G_id;
fprintf('FTLA inicial (sem compensador):\n');
disp(FTLA_int);

% Avaliar FTLA_int na frequência w_d
[mag_wd, phase_wd] = bode(FTLA_int, w_d);
mag_wd = squeeze(mag_wd);
phase_wd = squeeze(phase_wd);

fprintf('\nNa frequência w_d = %.2f rad/s:\n', w_d);
fprintf('|FTLA_int(w_d)| = %.4f\n', mag_wd);
fprintf('arg(FTLA_int(w_d)) = %.2f graus\n', phase_wd);

% Cálculo do zero do compensador (equação corrigida)
phase_required = M_p - 180 - phase_wd; % Fase necessária do compensador
if phase_required > 0
    phase_required = phase_required - 180;
end

fprintf('Fase necessária do compensador: %.2f graus\n', phase_required);

% Cálculo de w_z (zero do compensador)
if phase_required < -90
    w_z = w_d / tan(-phase_required * pi/180 - pi/2);
else
    w_z = w_d * tan(-phase_required * pi/180);
end

fprintf('Frequência do zero do compensador: w_z = %.2f rad/s\n', w_z);

% Cálculo do ganho k_d
k_d = 1 / (mag_wd * sqrt(1 + (w_d/w_z)^2));
fprintf('Ganho do compensador: k_d = %.4f\n', k_d);

% Controlador PI: C_i(s) = k_d * (s + w_z) / s
C_i = k_d * tf([1, w_z], [1, 0]);
fprintf('\nControlador PI C_i(s):\n');
disp(C_i);

% FTLA com compensador
FTLA_comp = FTLA_int * C_i;

% Malha fechada de corrente
G_cl_current = feedback(FTLA_comp, 1);
fprintf('Função de transferência de malha fechada de corrente:\n');
disp(G_cl_current);

%% 2. MALHA DE VELOCIDADE
fprintf('\n=== MALHA DE VELOCIDADE ===\n');

% Função de transferência da planta de velocidade
% G_n(w) = (30/pi) * (L_af * I_f) / (j * w * J)
% Convertendo para função de transferência no domínio s:
% G_n(s) = K / s, onde K = (30/pi) * (L_af * I_f) / J

K_vel = (30/pi) * (L_af * I_f) / J;
fprintf('Ganho da planta de velocidade: K_vel = %.4f\n', K_vel);

G_vel = tf(K_vel, [1, 0]);
fprintf('Função de transferência G_vel(s):\n');
disp(G_vel);

% Considerando k_s = 1 (ganho do sensor de velocidade)
k_s = 1;
FILA_vel = (1/k_s) * k_s * G_vel;  % FTLA para malha de velocidade
fprintf('FTLA para malha de velocidade:\n');
disp(FILA_vel);

%% 3. PLOTS DOS DIAGRAMAS DE BODE
fprintf('\n=== GERANDO DIAGRAMAS DE BODE ===\n');

% Frequência para análise
w = logspace(0, 6, 1000);  % 1 a 1e6 rad/s

% Figura 1: Malha de Corrente
figure('Position', [100, 100, 1200, 800]);

subplot(2,2,1);
[mag_int, phase_int] = bode(FTLA_int, w);
mag_int = squeeze(mag_int); phase_int = squeeze(phase_int);
semilogx(w, 20*log10(mag_int), 'b', 'LineWidth', 2); hold on;
[mag_comp, phase_comp] = bode(FTLA_comp, w);
mag_comp = squeeze(mag_comp); phase_comp = squeeze(phase_comp);
semilogx(w, 20*log10(mag_comp), 'r', 'LineWidth', 2);
grid on; title('Malha de Corrente - Magnitude'); 
xlabel('Frequência (rad/s)'); ylabel('Magnitude (dB)');
legend('FTLA_{inicial}', 'FTLA_{compensada}', 'Location', 'best');

subplot(2,2,2);
semilogx(w, phase_int, 'b', 'LineWidth', 2); hold on;
semilogx(w, phase_comp, 'r', 'LineWidth', 2);
grid on; title('Malha de Corrente - Fase');
xlabel('Frequência (rad/s)'); ylabel('Fase (graus)');
legend('FTLA_{inicial}', 'FTLA_{compensada}', 'Location', 'best');

subplot(2,2,3);
[mag_cl, phase_cl] = bode(G_cl_current, w);
mag_cl = squeeze(mag_cl); phase_cl = squeeze(phase_cl);
semilogx(w, 20*log10(mag_cl), 'g', 'LineWidth', 2);
grid on; title('Malha Fechada de Corrente - Magnitude');
xlabel('Frequência (rad/s)'); ylabel('Magnitude (dB)');

subplot(2,2,4);
semilogx(w, phase_cl, 'g', 'LineWidth', 2);
grid on; title('Malha Fechada de Corrente - Fase');
xlabel('Frequência (rad/s)'); ylabel('Fase (graus)');

% Figura 2: Malha de Velocidade
figure('Position', [100, 100, 800, 600]);

[mag_vel, phase_vel] = bode(FILA_vel, w);
mag_vel = squeeze(mag_vel); phase_vel = squeeze(phase_vel);

subplot(2,1,1);
semilogx(w, 20*log10(mag_vel), 'm', 'LineWidth', 2);
grid on; title('Malha de Velocidade - Magnitude');
xlabel('Frequência (rad/s)'); ylabel('Magnitude (dB)');

subplot(2,1,2);
semilogx(w, phase_vel, 'm', 'LineWidth', 2);
grid on; title('Malha de Velocidade - Fase');
xlabel('Frequência (rad/s)'); ylabel('Fase (graus)');

%% 4. ANÁLISE DE ESTABILIDADE
fprintf('\n=== ANÁLISE DE ESTABILIDADE ===\n');

% Margens de ganho e fase para malha de corrente compensada
[Gm_current, Pm_current, Wcg_current, Wcp_current] = margin(FTLA_comp);
fprintf('Malha de Corrente Compensada:\n');
fprintf('Margem de Ganho: %.2f dB em %.2f rad/s\n', 20*log10(Gm_current), Wcg_current);
fprintf('Margem de Fase: %.2f graus em %.2f rad/s\n', Pm_current, Wcp_current);

% Margens para malha de velocidade
[Gm_vel, Pm_vel, Wcg_vel, Wcp_vel] = margin(FILA_vel);
fprintf('\nMalha de Velocidade:\n');
fprintf('Margem de Ganho: %.2f dB em %.2f rad/s\n', 20*log10(Gm_vel), Wcg_vel);
fprintf('Margem de Fase: %.2f graus em %.2f rad/s\n', Pm_vel, Wcp_vel);

%% 5. PLOTS USANDO FUNÇÃO BODE DO MATLAB
figure('Position', [100, 100, 1000, 800]);

subplot(2,2,1);
bode(FTLA_int); grid on; title('FTLA Inicial de Corrente');
subplot(2,2,2);
bode(FTLA_comp); grid on; title('FTLA Compensada de Corrente');
subplot(2,2,3);
bode(G_cl_current); grid on; title('Malha Fechada de Corrente');
subplot(2,2,4);
bode(FILA_vel); grid on; title('FTLA Malha de Velocidade');

%% 6. RESPOSTA AO DEGRAU - MALHA FECHADA DE CORRENTE
figure('Position', [100, 100, 800, 400]);
t = 0:0.0001:0.02; % Tempo de 0 a 20ms
step(G_cl_current, t);
grid on; title('Resposta ao Degrau - Malha Fechada de Corrente');
xlabel('Tempo (s)'); ylabel('Amplitude');

% Informações da resposta ao degrau
stepinfo_current = stepinfo(G_cl_current);
fprintf('\nCaracterísticas da Resposta ao Degrau - Malha de Corrente:\n');
fprintf('Tempo de subida: %.4f s\n', stepinfo_current.RiseTime);
fprintf('Sobre-sinal: %.2f%%\n', stepinfo_current.Overshoot);
fprintf('Tempo de estabilização: %.4f s\n', stepinfo_current.SettlingTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
% closed loop voltage
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parâmetros do Sistema
V_in = 48;              % Tensão de barramento [V]
R_a = 0.016;            % Resistência de armadura [Ω]
L_a = 0.06;             % Indutância de armadura [H]
J = 0.018;              % Momento de inércia [kg·m²]
L_af = 0.078;           % Indutância mútua [H]
I_f = 1.6;              % Corrente de campo [A]
V_trip = 1;             % Valor de pico da triangular [V]

% Ganhos do controlador
k_i = 1;                % Ganho do sensor de corrente
kia = 1;                % Ganho do controlador de corrente
kpwm = 1/V_trip;        % Ganho do PWM

%% 1. MALHA DE CORRENTE (já projetada anteriormente)
fprintf('=== MALHA DE CORRENTE ===\n');

% Função de transferência da planta de corrente
G_id = tf(V_in, [L_a, R_a]);

% FTLA inicial (malha aberta sem compensador)
FTLA_int = k_i * kpwm * G_id;

% Parâmetros do controlador PI de corrente (do projeto anterior)
f_i = 1.66e4;           % Frequência do conversor [Hz]
f_u = f_i / 10;         % Frequência de cruzamento [Hz]
w_d = 2 * pi * f_u;     % Frequência angular de cruzamento [rad/s]
M_p = 60 * pi/180;      % 60 graus em radianos

% Avaliar FTLA_int na frequência w_d
[mag_wd, phase_wd] = bode(FTLA_int, w_d);
mag_wd = squeeze(mag_wd);
phase_wd = squeeze(phase_wd);

% Cálculo do controlador PI de corrente
phase_required = M_p - 180 - phase_wd;
w_z = w_d / tan(-phase_required * pi/180 - pi/2);
k_d = 1 / (mag_wd * sqrt(1 + (w_d/w_z)^2));

% Controlador PI de corrente
C_i = k_d * tf([1, w_z], [1, 0]);

% FTLA com compensador de corrente
FTLA_comp = FTLA_int * C_i;

% Malha fechada de corrente
G_cl_current = feedback(FTLA_comp, 1);

%% 2. MALHA DE TENSÃO - Projeto do Controlador PI
fprintf('\n=== PROJETO DO CONTROLADOR PI PARA MALHA DE TENSÃO ===\n');

% Frequências de projeto para tensão
f_d = 1.66e3;           % Frequência de cruzamento da corrente [Hz]
f_0 = f_d / 100;        % Frequência de cruzamento da tensão [Hz]
w_c1 = 2 * pi * f_0;    % Frequência angular de cruzamento [rad/s]

fprintf('Frequência de cruzamento da tensão: f_0 = %.2f Hz\n', f_0);
fprintf('Frequência angular de cruzamento: w_c1 = %.2f rad/s\n', w_c1);

% Margem de fase desejada
M_b = 60 * pi/180;      % 60 graus em radianos
fprintf('Margem de fase desejada: M_b = %.2f graus\n', M_b*180/pi);

% 2.1 Função de transferência da planta de tensão
% Considerando que a tensão é controlada através da corrente de campo
% Para um motor CC: V_t = K_v * w_m + R_a * I_a
% Onde K_v = L_af * I_f

K_v = L_af * I_f;       % Constante de tensão
fprintf('Constante de tensão K_v = L_af * I_f = %.4f\n', K_v);

% Supondo uma dinâmica de primeira ordem para a tensão
% Com constante de tempo τ_v (a ser determinada)
tau_v = 0.01;           % Constante de tempo estimada [s]
G_v = tf(K_v, [tau_v, 1]);  % Planta de tensão

fprintf('Função de transferência da planta de tensão G_v(s):\n');
disp(G_v);

% 2.2 FTLA inicial para tensão (considerando malha de corrente fechada)
% A malha de tensão vê a malha de corrente como um bloco
FTLA_v_int = G_cl_current * G_v;
fprintf('FTLA inicial para tensão (com malha de corrente fechada):\n');
disp(FTLA_v_int);

% 2.3 Avaliar FTLA_v_int na frequência w_c1
[mag_wc1, phase_wc1] = bode(FTLA_v_int, w_c1);
mag_wc1 = squeeze(mag_wc1);
phase_wc1 = squeeze(phase_wc1);

fprintf('\nNa frequência w_c1 = %.2f rad/s:\n', w_c1);
fprintf('|FTLA_v_int(w_c1)| = %.4f\n', mag_wc1);
fprintf('arg(FTLA_v_int(w_c1)) = %.2f graus\n', phase_wc1);

% 2.4 Cálculo do compensador PI para tensão
% avgPIL = (M_b - pi - arg(FTLA_v_int(w_c1))) * 180/pi
avgPIL = (M_b - pi - phase_wc1*pi/180) * 180/pi;
fprintf('avgPIL = %.2f graus\n', avgPIL);

% Cálculo do polo do compensador (ω_p)
phase_comp_req = M_b - 180 - phase_wc1; % Fase necessária do compensador
fprintf('Fase necessária do compensador: %.2f graus\n', phase_comp_req);

if phase_comp_req > -90
    w_p = w_c1 * tan(phase_comp_req * pi/180);
else
    w_p = w_c1 / tan(-phase_comp_req * pi/180 - pi/2);
end

fprintf('Frequência do polo do compensador: w_p = %.4f rad/s\n', w_p);

% Cálculo do ganho k_o
k_o = 1 / (mag_wc1 * sqrt(1 + (w_c1/w_p)^2));
fprintf('Ganho do compensador: k_o = %.6f\n', k_o);

% 2.5 Controlador PI para tensão: C_v(s) = k_o * (s + w_p) / s
C_v = k_o * tf([1, w_p], [1, 0]);
fprintf('\nControlador PI para tensão C_v(s):\n');
disp(C_v);

% 2.6 FTLA com compensador de tensão
FTLA_v_comp = FTLA_v_int * C_v;

% 2.7 Malha fechada de tensão
G_cl_voltage = feedback(FTLA_v_comp, 1);
fprintf('Função de transferência de malha fechada de tensão:\n');
disp(G_cl_voltage);

%% 3. PLOTS DOS DIAGRAMAS DE BODE - MALHA DE TENSÃO
fprintf('\n=== DIAGRAMAS DE BODE - MALHA DE TENSÃO ===\n');

% Frequência para análise
w = logspace(-1, 4, 1000);  % 0.1 a 10000 rad/s

figure('Position', [100, 100, 1200, 800]);

% Subplot 1: Magnitude
subplot(2,2,1);
[mag_v_int, phase_v_int] = bode(FTLA_v_int, w);
mag_v_int = squeeze(mag_v_int); phase_v_int = squeeze(phase_v_int);
[mag_v_comp, phase_v_comp] = bode(FTLA_v_comp, w);
mag_v_comp = squeeze(mag_v_comp); phase_v_comp = squeeze(phase_v_comp);

semilogx(w, 20*log10(mag_v_int), 'b', 'LineWidth', 2); hold on;
semilogx(w, 20*log10(mag_v_comp), 'r', 'LineWidth', 2);
grid on; title('Malha de Tensão - Magnitude'); 
xlabel('Frequência (rad/s)'); ylabel('Magnitude (dB)');
legend('FTLA_{v inicial}', 'FTLA_{v compensada}', 'Location', 'best');

% Subplot 2: Fase
subplot(2,2,2);
semilogx(w, phase_v_int, 'b', 'LineWidth', 2); hold on;
semilogx(w, phase_v_comp, 'r', 'LineWidth', 2);
grid on; title('Malha de Tensão - Fase');
xlabel('Frequência (rad/s)'); ylabel('Fase (graus)');
legend('FTLA_{v inicial}', 'FTLA_{v compensada}', 'Location', 'best');

% Subplot 3: Malha Fechada de Tensão - Magnitude
subplot(2,2,3);
[mag_cl_v, phase_cl_v] = bode(G_cl_voltage, w);
mag_cl_v = squeeze(mag_cl_v); phase_cl_v = squeeze(phase_cl_v);
semilogx(w, 20*log10(mag_cl_v), 'g', 'LineWidth', 2);
grid on; title('Malha Fechada de Tensão - Magnitude');
xlabel('Frequência (rad/s)'); ylabel('Magnitude (dB)');

% Subplot 4: Malha Fechada de Tensão - Fase
subplot(2,2,4);
semilogx(w, phase_cl_v, 'g', 'LineWidth', 2);
grid on; title('Malha Fechada de Tensão - Fase');
xlabel('Frequência (rad/s)'); ylabel('Fase (graus)');

%% 4. PLOTS COMPARATIVOS - TODAS AS MALHAS
figure('Position', [100, 100, 1200, 600]);

% Magnitude de todas as malhas
subplot(1,2,1);
[mag_cl_curr, ~] = bode(G_cl_current, w);
mag_cl_curr = squeeze(mag_cl_curr);
semilogx(w, 20*log10(mag_cl_curr), 'b', 'LineWidth', 2); hold on;
semilogx(w, 20*log10(mag_cl_v), 'r', 'LineWidth', 2);
grid on; title('Comparação - Magnitude das Malhas Fechadas');
xlabel('Frequência (rad/s)'); ylabel('Magnitude (dB)');
legend('Malha de Corrente', 'Malha de Tensão', 'Location', 'best');

% Fase de todas as malhas
subplot(1,2,2);
[~, phase_cl_curr] = bode(G_cl_current, w);
phase_cl_curr = squeeze(phase_cl_curr);
semilogx(w, phase_cl_curr, 'b', 'LineWidth', 2); hold on;
semilogx(w, phase_cl_v, 'r', 'LineWidth', 2);
grid on; title('Comparação - Fase das Malhas Fechadas');
xlabel('Frequência (rad/s)'); ylabel('Fase (graus)');
legend('Malha de Corrente', 'Malha de Tensão', 'Location', 'best');

%% 5. ANÁLISE DE ESTABILIDADE - MALHA DE TENSÃO
fprintf('\n=== ANÁLISE DE ESTABILIDADE - MALHA DE TENSÃO ===\n');

% Margens de ganho e fase para malha de tensão compensada
[Gm_voltage, Pm_voltage, Wcg_voltage, Wcp_voltage] = margin(FTLA_v_comp);
fprintf('Malha de Tensão Compensada:\n');
fprintf('Margem de Ganho: %.2f dB em %.2f rad/s\n', 20*log10(Gm_voltage), Wcg_voltage);
fprintf('Margem de Fase: %.2f graus em %.2f rad/s\n', Pm_voltage, Wcp_voltage);

% Verificação na frequência de cruzamento w_c1
[mag_check, phase_check] = bode(FTLA_v_comp, w_c1);
mag_check = squeeze(mag_check);
fprintf('\nVerificação em w_c1 = %.2f rad/s:\n', w_c1);
fprintf('|FTLA_v_comp(w_c1)| = %.4f (0 dB esperado)\n', mag_check);
fprintf('Ganho em dB: %.2f dB\n', 20*log10(mag_check));

%% 6. RESPOSTA AO DEGRAU - MALHA FECHADA DE TENSÃO
figure('Position', [100, 100, 800, 400]);
t = 0:0.001:0.1; % Tempo de 0 a 100ms
step(G_cl_voltage, t);
grid on; title('Resposta ao Degrau - Malha Fechada de Tensão');
xlabel('Tempo (s)'); ylabel('Amplitude');

% Informações da resposta ao degrau
stepinfo_voltage = stepinfo(G_cl_voltage);
fprintf('\nCaracterísticas da Resposta ao Degrau - Malha de Tensão:\n');
fprintf('Tempo de subida: %.4f s\n', stepinfo_voltage.RiseTime);
fprintf('Sobre-sinal: %.2f%%\n', stepinfo_voltage.Overshoot);
fprintf('Tempo de estabilização: %.4f s\n', stepinfo_voltage.SettlingTime);

%% 7. PLOTS USANDO FUNÇÃO BODE DO MATLAB
figure('Position', [100, 100, 1000, 800]);

subplot(2,2,1);
bode(FTLA_v_int); grid on; title('FTLA Inicial de Tensão');
subplot(2,2,2);
bode(FTLA_v_comp); grid on; title('FTLA Compensada de Tensão');
subplot(2,2,3);
bode(G_cl_voltage); grid on; title('Malha Fechada de Tensão');
subplot(2,2,4);
margin(FTLA_v_comp); grid on; title('Margens de Estabilidade - Tensão');

fprintf('\n=== PROJETO CONCLUÍDO ===\n');
fprintf('Controlador PI de Tensão projetado com sucesso!\n');
fprintf('Frequência de cruzamento: %.2f Hz\n', f_0);
fprintf('Margem de fase obtida: %.2f graus\n', Pm_voltage);

