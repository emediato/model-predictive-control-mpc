%% MATLAB Code for Induction Motor Control System Specifications
% TRABALHO 03 - mit
% Sistema de Controle Vetorial com Inversor Trifásico

clear all; clc; close all;

%% DEFINIÇÃO DE UNIDADES E CONSTANTES
k = 1e3;    % Kilo (10^3)
m = 1e-3;   % Mili (10^-3)
j = sqrt(-1); % un imaginária

% Constantes de conversão
k_rads_prpm = 2 * pi / 60; % Conversão RPM para rad/s

%%  ESPECIFICAÇÕES DO MOTOR DE INDUÇÃO
% Motor: 3hp = 2.2371kW
fprintf(' ESPECIFICAÇÕES DO MOTOR \n');

% Potência Nominal
P_nom = 2.2371 * k; % 2.2371 kW = 2237.1 W
fprintf('Potência Nominal: P_nom = %.4f W\n', P_nom);

% Parâmetros do Motor
B = 0;              % Coeficiente de atrito [N·m·s/rad]
I_d_nom = 8;        % Corrente direta nominal [A]
Ls = 2 * m;         % Indutância do estator [H] (2 mH)
Lr = 2 * m;         % Indutância do rotor [H] (2 mH)
Lm = 69.3 * m;      % Indutância de magnetização [H] (69.3 mH)
J = 89 * m;         % Momento de inércia [kg·m²] (89 mkg)
P = 4;              % Número de polos
Rs = 0.435;         % Resistência do estator [Ω]
Rr = 0.816;         % Resistência do rotor [Ω]

% Exibir parâmetros do motor
motor_params = table([P; Rs; Rr; Ls*1000; Lr*1000; Lm*1000; J*1000; B; I_d_nom], ...
    'VariableNames', {'Valor'}, ...
    'RowNames', {'Número de polos', 'Rs [Ω]', 'Rr [Ω]', 'Ls [mH]', 'Lr [mH]', ...
                 'Lm [mH]', 'J [mkg]', 'B [N·m·s/rad]', 'I_d_nom [A]'});
disp(motor_params);

%% CÁLCULO DOS COEFICIENTES DO MOTOR
fprintf('\n COEFICIENTES CALCULADOS \n');

% Cálculo das indutâncias equivalentes
L_ss = Ls + Lm;     % Indutância total do estator
L_rr = Lr + Lm;     % Indutância total do rotor
L_a = L_ss - (Lm^2 / L_rr); % Indutância de eixo
R_a = Rs + (L_ss / L_rr) * Rr; % Resistência equivalente
T_a = L_a / R_a;    % Constante de tempo elétrica
K_a = 1 / R_a;      % Ganho estático

% Constante de torque
K_t = (3/2) * (P/2) * (Lm^2 / L_rr) * I_d_nom;

% Exibir coeficientes calculados
coefficients = table([L_ss; L_rr; L_a*1000; R_a; K_a; T_a*1000; K_t], ...
    'VariableNames', {'Valor'}, ...
    'RowNames', {'L_ss [H]', 'L_rr [H]', 'L_a [mH]', 'R_a [Ω]', ...
                 'K_a [1/Ω]', 'T_a [ms]', 'K_t [N·m/A]'});
disp(coefficients);

%% ESPECIFICAÇÕES DO SISTEMA DE CONVERSÃO
fprintf('\n ESPECIFICAÇÕES DO SISTEMA \n');

% Tensão DC-link
V_dc = 500;         % Tensão da bateria/DC-link [V]

% Frequência de comutação
f_s = 40 * k;       % 40 kHz = 40000 Hz

% Valor de pico da portadora triangular
V_tri = 5;          % [V]

% Ganhos do sistema
k_r = V_dc / 2;     % Ganho do inversor [V]
k_tig = 1;          % Ganho do sensor de corrente direta
k_td = 1;           % Ganho do sensor de corrente quadratura
k_v = 1;            % Ganho do sensor de tensão
k_pum = 1 / V_tri;  % Ganho do modulador PWM

% Atraso de comutação
tau_r = (1/2) * (1/f_s); % [s]

% Função de transferência do inversor
s = tf('s');
G_inv = k_r / (1 + tau_r * s); % Modelo de atraso de 1ª ordem
G_inv_pade = k_r * (1 - tau_r*s/2) / (1 + tau_r*s/2); % Aproximação de Padé 
G_inv_delay = k_r * exp(-tau_r * s); % Atraso puro (para simulação no domínio do tempo)

% Exibir especificações do sistema
system_specs = table([V_dc; f_s/1000; V_tri; k_r; k_tig; k_td; k_v; k_pum; tau_r*1e6], ...
    'VariableNames', {'Valor'}, ...
    'RowNames', {'V_dc [V]', 'f_s [kHz]', 'V_tri [V]', 'k_r [V]', ...
                 'k_tig', 'k_td', 'k_v', 'k_pum', 'tau_r [μs]'});
disp(system_specs);

%% ANÁLISE EM FREQUÊNCIA
fprintf('\n ANÁLISE EM FREQUÊNCIA \n');

% Definir vetor de frequências
f = logspace(0, 5, 1000); % De 1 Hz a 100 kHz
w = 2 * pi * f;           % Frequência angular [rad/s]

% Resposta em frequência do inversor
[mag_inv, phase_inv] = bode(G_inv, w);
mag_inv = squeeze(mag_inv);
phase_inv = squeeze(phase_inv);
% 
% % Plotar resposta em frequência
% figure('Position', [100, 100, 800, 600]);
% 
% subplot(2,1,1);
% semilogx(f, 20*log10(mag_inv), 'b', 'LineWidth', 2);
% grid on;
% xlabel('Frequência [Hz]', 'FontSize', 12);
% ylabel('Magnitude [dB]', 'FontSize', 12);
% title('Resposta em Frequência do Inversor ', 'FontSize', 14);
% xlim([1, 1e5]);
% ylim([20*log10(k_r)-10, 20*log10(k_r)+5]);
% 
% subplot(2,1,2);
% semilogx(f, phase_inv, 'r', 'LineWidth', 2);
% grid on;
% xlabel('Frequência [Hz]', 'FontSize', 12);
% ylabel('Fase [graus]', 'FontSize', 12);
% title('Resposta de Fase do Inversor', 'FontSize', 14);
% xlim([1, 1e5]);


%%******************************************************************************
%%******************************************************************************

%%******************************************************************************
%%******************************************************************************

%%******************************************************************************
%%******************************************************************************
%%  FUNÇÃO DE TRANSFERÊNCIA DA PLANTA (EIXO q)
fprintf('\n FUNÇÃO DE TRANSFERÊNCIA DA PLANTA \n');


% Função de transferência da planta G_iq(s)
s = tf('s');
G_iq = K_a  / (1 + s * T_a)  % Planta de primeira ordem

% Simplificação: polos e zeros
[p_iq, z_iq] = pzmap(G_iq);
fprintf('Polo da planta: s = %.2f rad/s\n', p_iq);
fprintf('Equivalente: s = %.2f rad/s \n', 1/T_a);

num_iq = K_a  ;
den_iq = [T_a, 1];
G_iq_alt = tf(num_iq, den_iq);

% Verificação
fprintf('\nVerificação da planta:\n');
fprintf('K_a/(1+sT_a) = %.3f/(1 + s*%.3f)\n', K_a, T_a);
fprintf('Polo: s = -1/T_a = -%.2f rad/s\n', 1/T_a);

%% FUNÇÃO DE TRANSFERÊNCIA DE LAÇO ABERTO NÃO COMPENSADO
fprintf('\n LAÇO ABERTO NÃO COMPENSADO \n');
k_iq = 1;           % Ganho do sensor de corrente (eixo q)

% FTLA não compensada (sem controlador)
FTLA_NC_iq = k_pum * k_r * G_iq * k_iq;

% Simplificar para a forma da imagem
FTLA_NC_iq_simplified = tf(k_pum * k_r * K_a * k_iq, [T_a, 1]);

% Mostrar função de transferência
fprintf('FTLA_NC_iq(s) = k_pum * k_r * G_iq(s) * k_iq\n');
fprintf('FTLA_NC_iq(s) = %.3f * %.0f * %.3f/(1+s*%.3f) * %.0f\n', ...
    k_pum, k_r, K_a, T_a, k_iq);
fprintf('FTLA_NC_iq(s) = %.4f/(s + %.2f)\n', ...
    k_pum * k_r * K_a * k_iq / T_a, 1/T_a);


%% ESPECIFICAÇÕES DO CONTROLADOR
fprintf('\n ESPECIFICAÇÕES DO CONTROLADOR PI \n');

% Frequência de cruzamento desejada
f_c_iq = f_s / 10;  % f_s/10 = 4 kHz
w_c_iq = 2 * pi * f_c_iq;  % Frequência angular de cruzamento

% Margem de fase desejada
MF_iq_deg = 60;     % Margem de fase em graus
MF_iq_rad = deg2rad(MF_iq_deg);  % Converter para radianos

fprintf('Frequência de cruzamento: f_c_iq = %.0f Hz\n', f_c_iq);
fprintf('                        w_c_iq = %.3e rad/s\n', w_c_iq);
fprintf('Margem de fase desejada: MF = %.0f° (%.3f rad)\n', MF_iq_deg, MF_iq_rad);

%ANÁLISE DO SISTEMA NÃO COMPENSADO NA FREQUÊNCIA DE CRUZAMENTO
fprintf('\n ANÁLISE NA FREQUÊNCIA DE CRUZAMENTO \n');

% Avaliar FTLA_NC na frequência de cruzamento
[mag_NC, phase_NC] = bode(FTLA_NC_iq, w_c_iq);
mag_NC_db = 20*log10(mag_NC);  % Magnitude em dB
phase_NC_deg = phase_NC;       % Fase em graus

mag_nc_at_wc = mag_NC;

% Resposta em frequência da FTLA_NC_IQ
[mag_FTLA_NC_iq, phase_FTLA_NC_iq] = bode(FTLA_NC_iq, w_c_iq);


phase_nc_at_wc_deg = phase_FTLA_NC_iq;
phase_nc_at_wc_rad = deg2rad(phase_FTLA_NC_iq)

mag_FTLA_NC_iq = squeeze(mag_FTLA_NC_iq);
phase_FTLA_NC_iq = squeeze(phase_FTLA_NC_iq);
% 
% fprintf('a) Avaliação de FTLA_NC_id(jω_c_id):\n');
% fprintf('   |FTLA_NC_iq(jω_c_iq)| = %.8f\n', mag_FTLA_NC_iq);
% fprintf('   ∠FTLA_NC_iq(jω_c_iq) = %.6f°\n', phase_nc_at_wc_deg);
% fprintf('   ∠FTLA_NC_iq(jω_c_iq) = %.6f rad\n', phase_nc_at_wc_rad);

%% CONTROLADOR PI
fprintf('\nb) Cálculo do ângulo necessário do compensador:\n');
fprintf('   ângulo_PI = MF_id_rad - pi - ∠FTLA_NC_id(jω_c_id)\n');
fprintf('   ângulo_PI = %.6f - %.6f - %.6f\n', ...
    MF_iq_rad, pi, phase_nc_at_wc_rad);
ang_PI_rad = MF_iq_rad - pi - phase_nc_at_wc_rad;
ang_PI_deg = rad2deg(ang_PI_rad);
fprintf('   ângulo_PI = %.6f rad = %.3f°\n', ang_PI_rad, ang_PI_deg);


% Converter ângulo_PI para graus para usar função tand()
ang_PI_deg_for_tan = ang_PI_deg;

% % Calcular ω_z usando a tangente
% fprintf('\nb) Cálculo de ω_z:\n');
% fprintf('   ω_z = ω_c_id / tan(ângulo_PI + 90°)\n');
% fprintf('   ω_z = %.3e / tan(%.3f° + 90°)\n', w_c_iq, ang_PI_deg_for_tan);
% fprintf('   ω_z = %.3e / tan(%.3f°)\n', w_c_iq, ang_PI_deg_for_tan + 90);

% ω_z = ω_c / tan(MF_rad - pi/2 - ∠FTLA_NC(ω_c))
arg_for_tan = MF_iq_rad - pi/2 - deg2rad(phase_NC);
%arg_for_tan = MF_iq_rad - pi/2 - phase_nc_at_wc_rad;

w_z_iq = w_c_iq / tan(arg_for_tan);
tau_z_iq = 1/w_z_iq;

% fprintf('   ω_z = ω_c_id / tan(MF_id_rad - pi/2 - ∠FTLA_NC(ω_c_id))\n');
% fprintf('   arg = %.6f - %.6f - %.6f\n', MF_iq_rad, pi/2, phase_nc_at_wc_rad);
% fprintf('   arg = %.6f rad\n', arg_for_tan);
% fprintf('   ω_z = %.3e / tan(%.6f)\n', w_c_iq, arg_for_tan);
fprintf('   ω_z = %.6e rad/s\n', w_z_iq);

mag_pi_at_wc_expr = sqrt(w_c_iq^2 + w_z_iq^2) / w_c_iq;
fprintf('   |C(jω_c)| = Kc × %.6f\n', mag_pi_at_wc_expr);


% Condição de módulo completa
fprintf('\nc) Resolvendo para Kc:\n');
fprintf('   Kc × %.6f × %.8f = 1\n', mag_pi_at_wc_expr, mag_nc_at_wc);
fprintf('   Kc = 1 / (%.6f × %.8f)\n', mag_pi_at_wc_expr, mag_nc_at_wc);

% Calcular Kc
Kc_iq = 1 / (mag_pi_at_wc_expr * mag_nc_at_wc);
fprintf('   Kc = %.8f\n', Kc_iq);



% Compensador PI: C(s) = Kc × (s + ω_z)/s
C_pi_iq = Kc_iq * (s + w_z_iq) / s;

% Forma alternativa: C(s) = Kp + Ki/s
Kp_iq = Kc_iq;
Ki_iq = Kc_iq * w_z_iq;


G_iq = K_a / (1 + s * T_a);

% FTLA não compensada
FTLA_NC_iq = k_pum * k_r * G_iq * k_iq;

% FTLA compensada
FTLA_C_iq = C_pi_iq * FTLA_NC_iq


%%

% % Plotar resposta em frequência
% figure('Position', [100, 100, 800, 600]);
% 
% subplot(2,1,1);
% semilogx(f, 20*log10(mag_FTLA_NC_iq), 'b', 'LineWidth', 2);
% grid on;
% xlabel('Frequência [Hz]', 'FontSize', 12);
% ylabel('Magnitude [dB]', 'FontSize', 12);
% title('Resposta em Frequência do Inversor FTLA_NC_iq', 'FontSize', 14);
% % xlim([1, 1e5]);
% % ylim([20*log10(k_r)-10, 20*log10(k_r)+5]);
% 
% subplot(2,1,2);
% semilogx(f, phase_FTLA_NC_iq, 'r', 'LineWidth', 2);
% grid on;
% xlabel('Frequência [Hz]', 'FontSize', 12);
% ylabel('Fase [graus]', 'FontSize', 12);
% title('Resposta de Fase do Inversor', 'FontSize', 14);
% 

% FTLA não compensada
[mag_NC, phase_NC] = bode(FTLA_NC_iq, w);
mag_NC = squeeze(mag_NC);
phase_NC = squeeze(phase_NC);
mag_NC_dB = 20 * log10(mag_NC);     % Magnitude em dB

% FTLA compensada
[mag_C, phase_C] = bode(FTLA_C_iq, w);
mag_C = squeeze(mag_C);
phase_C = squeeze(phase_C);
mag_C_dB = 20 * log10(mag_C);       % Magnitude em dB

    
% Calcular frequência de cruzamento real (onde magnitude = 0 dB)
[~, iqx_cross] = min(abs(mag_C_dB));
f_c_real = f(iqx_cross);

% Calcular margem de fase na frequência de cruzamento
[~, phase_at_fc] = bode(FTLA_C_iq, 2*pi*f_c_real);
phase_at_fc = squeeze(phase_at_fc);
MF_real = 180 + phase_at_fc;  % Margem de fase real



figure('Position', [50, 50, 1400, 600], 'Name', 'Diagramas de Bode Giq - FTLA Não Compensada vs Compensada');
% Subplot 1: Magnitude
subplot(1,2,1);
semilogx(f, mag_NC_dB, 'b-', 'LineWidth', 2, 'DisplayName', 'FTLA_{NC} IQ - Não Compensada');
hold on;
semilogx(f, mag_C_dB, 'r-', 'LineWidth', 2, 'DisplayName', 'FTLA_{C}  IQ- Compensada');
grid on;

% Linhas de referência
yline(0, 'k--', 'LineWidth', 1, 'DisplayName', '0 dB');
if exist('f_c', 'var')
    xline(f_c, 'g--', 'LineWidth', 1.5, 'DisplayName', sprintf('f_c = %.0f Hz', f_c));
end
xline(f_c_real, 'm--', 'LineWidth', 1.5, 'DisplayName', sprintf('f_c real = %.0f Hz', f_c_real));

xlabel('Frequência [Hz]', 'FontSize', 12);
ylabel('Magnitude [dB]', 'FontSize', 12);
title('Comparação de Magnitude  IQ ', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);

ylim([-60, 60]);


% Subplot 2: Fase
subplot(1,2,2);
semilogx(f, phase_NC, 'b-', 'LineWidth', 2, 'DisplayName', 'FTLA_{NC} IQ - Não Compensada');
hold on;
semilogx(f, phase_C, 'r-', 'LineWidth', 2, 'DisplayName', 'FTLA_{C} IQ - Compensada');
grid on;

% Linhas de referência
yline(-180, 'k--', 'LineWidth', 1, 'DisplayName', '-180°');
yline(-180+60, 'g--', 'LineWidth', 1.5, 'DisplayName', 'MF = 60°');
if exist('f_c', 'var')
    xline(f_c, 'g--', 'LineWidth', 1.5);
end
xline(f_c_real, 'm--', 'LineWidth', 1.5);

xlabel('Frequência [Hz]', 'FontSize', 12);
ylabel('Fase [graus]', 'FontSize', 12);
title('Comparação de Fase IQ', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);

ylim([-270, 0])



%%******************************************************************************
%%******************************************************************************

%%******************************************************************************
%%******************************************************************************

%%******************************************************************************
%%******************************************************************************

%%  FUNÇÃO DE TRANSFERÊNCIA DA PLANTA CORRENTE DIRETA (EIXO d)
G_id = G_iq 



%%******************************************************************************
%%******************************************************************************

%%******************************************************************************
%%******************************************************************************

%%******************************************************************************
%%******************************************************************************
%%  FUNÇÃO DE TRANSFERÊNCIA DA PLANTA VELOCIDADE - MALHA EXTERNA

% Especificações da malha de velocidade
MF_mi_deg = 60;         % Margem de fase desejada: 60 graus
f_c_mi = f_c_iq / 10;   % Frequência de cruzamento: f_c_iq/10 = 400 Hz

%Frequência angular de cruzamento
w_c_mi = 2 * pi * f_c_mi

% 2.2 Margem de fase em radianos
MF_mi_rad = deg2rad(MF_mi_deg)


s = tf('s');
K_i = K_t;  % Constante de conversão corrente-torque

% Função de transferência da velocidade (rad/s) por corrente (A)
G_mi_rad = K_i / (J*s + B);

% Para converter de rad/s para RPM:  rad/s × (60/(2π)) = RPM
conversion_factor = 60/(2*pi);  % 30/π ≈ 9.5493

% Função de transferência da velocidade (RPM) por corrente (A)
G_mi = G_mi_rad * conversion_factor;

% FTLA_NC_mi = k_v * G_mi(s)
FTLA_NC_mi = k_v * G_mi;

% Mostrar a função
[num_nc_mi, den_nc_mi] = tfdata(FTLA_NC_mi, 'v');
fprintf('FTLA_NC_mi(s) = k_v × G_mi(s)\n');
fprintf('FTLA_NC_mi(s) = %.0f × G_mi(s)\n', k_v);

[mag_nc_at_wc_mi, phase_nc_at_wc_mi] = bode(FTLA_NC_mi, w_c_mi);
phase_nc_at_wc_mi_deg = phase_nc_at_wc_mi;
phase_nc_at_wc_mi_rad = deg2rad(phase_nc_at_wc_mi_deg);

fprintf('   |FTLA_NC_mi(jω_c_mi)| = %.8f\n', mag_nc_at_wc_mi);
fprintf('   ∠FTLA_NC_mi(jω_c_mi) = %.6f°\n', phase_nc_at_wc_mi_deg);
fprintf('   ∠FTLA_NC_mi(jω_c_mi) = %.6f rad\n', phase_nc_at_wc_mi_rad);


% Para um integrador puro (K/s), a fase é sempre -90°
expected_phase = -90;

ang_PI_mi_rad = MF_mi_rad - pi - phase_nc_at_wc_mi_rad;
ang_PI_mi_deg = rad2deg(ang_PI_mi_rad)


arg_for_tan_mi = MF_mi_rad - pi/2 - phase_nc_at_wc_mi_rad;
w_z_mi = w_c_mi / tan(arg_for_tan_mi);

fprintf('   ω_z_mi = ω_c_mi / tan(arg)\n');
fprintf('   ω_z_mi = %.3e / tan(%.6f)\n', w_c_mi, arg_for_tan_mi);
fprintf('   ω_z_mi = %.6e rad/s\n', w_z_mi);


part1 = w_c_mi / sqrt(w_c_mi^2 + w_z_mi^2);
part2 = 1 / mag_nc_at_wc_mi;
Kc_mi = part1 * part2

tau_z_mi = 1 / w_z_mi;


% Compensador PI: C(s) = Kc × (s + ω_z)/s
C_pi_mi = Kc_mi * (s + w_z_mi) / s;

% Forma alternativa: C(s) = Kp + Ki/s
Kp_mi = Kc_mi;
Ki_mi = Kc_mi * w_z_mi;


FTLA_C_mi = C_pi_mi * FTLA_NC_mi

fprintf('FTLA_C_mi(s) = C_pi_mi(s) × FTLA_NC_mi(s)\n');


%%
% Simplificar a expressão

% Definir vetor de frequências (0.1 Hz a 10 kHz)
f_min_mi = 0.1;        % 0.1 Hz
f_max_mi = 10e3;       % 10 kHz
N_points_mi = 1000;

f_mi = logspace(log10(f_min_mi), log10(f_max_mi), N_points_mi);
w_mi = 2 * pi * f_mi;

% 12.1 Calcular respostas em frequência
fprintf('Calculando respostas em frequência...\n');

% FTLA não compensada
[mag_NC_mi, phase_NC_mi] = bode(FTLA_NC_mi, w_mi);
mag_NC_mi = squeeze(mag_NC_mi);
phase_NC_mi = squeeze(phase_NC_mi);
mag_NC_mi_dB = 20 * log10(mag_NC_mi);

% FTLA compensada
[mag_C_mi, phase_C_mi] = bode(FTLA_C_mi, w_mi);
mag_C_mi = squeeze(mag_C_mi);
phase_C_mi = squeeze(phase_C_mi);
mag_C_mi_dB = 20 * log10(mag_C_mi);

% 12.2 Calcular frequência de cruzamento real
[~, idx_cross_mi] = min(abs(mag_C_mi_dB));
f_c_mi_real = f_mi(idx_cross_mi);

% 12.3 Calcular margem de fase real
[~, phase_at_fc_mi] = bode(FTLA_C_mi, 2*pi*f_c_mi_real);
phase_at_fc_mi = squeeze(phase_at_fc_mi);
MF_mi_real = 180 + phase_at_fc_mi;

%COMPARAÇÃO FTLA NC vs C (VELOCIDADE)
figure('Position', [50, 50, 1400, 600], 'Name', 'Malha de Velocidade - Bode');

% Subplot 1: Magnitude
subplot(1,2,1);
semilogx(f_mi, mag_NC_mi_dB, 'b-', 'LineWidth', 2, 'DisplayName', 'FTLA_{NC,mi} - Não Compensada');
hold on;
semilogx(f_mi, mag_C_mi_dB, 'r-', 'LineWidth', 2, 'DisplayName', 'FTLA_{C,mi} - Compensada');
grid on;

% Linhas de referência
yline(0, 'k--', 'LineWidth', 1, 'DisplayName', '0 dB');
xline(f_c_mi, 'g--', 'LineWidth', 1.5, 'DisplayName', sprintf('f_c desejada = %.0f Hz', f_c_mi));
xline(f_c_mi_real, 'm--', 'LineWidth', 1.5, 'DisplayName', sprintf('f_c real = %.1f Hz', f_c_mi_real));

xlabel('Frequência [Hz]', 'FontSize', 12);
ylabel('Magnitude [dB]', 'FontSize', 12);
title('Malha de Velocidade - Magnitude', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
xlim([f_min_mi, f_max_mi]);
ylim([-60, 60]);

% Subplot 2: Fase
subplot(1,2,2);
semilogx(f_mi, phase_NC_mi, 'b-', 'LineWidth', 2, 'DisplayName', 'FTLA_{NC,mi} - Não Compensada');
hold on;
semilogx(f_mi, phase_C_mi, 'r-', 'LineWidth', 2, 'DisplayName', 'FTLA_{C,mi} - Compensada');
grid on;

% Linhas de referência
yline(-180, 'k--', 'LineWidth', 1, 'DisplayName', '-180°');
yline(-180+MF_mi_deg, 'g--', 'LineWidth', 1.5, 'DisplayName', sprintf('MF = %.0f°', MF_mi_deg));
xline(f_c_mi, 'g--', 'LineWidth', 1.5);
xline(f_c_mi_real, 'm--', 'LineWidth', 1.5);

xlabel('Frequência [Hz]', 'FontSize', 12);
ylabel('Fase [graus]', 'FontSize', 12);
title('Malha de Velocidade - Fase', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
xlim([f_min_mi, f_max_mi]);
ylim([-270, 0]);



%% CÁLCULO DE PARÂMETROS ADICIONAIS
fprintf('\n PARÂMETROS ADICIONAIS \n');

% Velocidade síncrona nominal (considerando 60 Hz)
f_e = 60; % Frequência elétrica nominal [Hz]
w_e = 2 * pi * f_e; % Velocidade elétrica [rad/s]
w_sync = w_e / (P/2); % Velocidade síncrona [rad/s]
N_sync = w_sync / k_rads_prpm; % Velocidade síncrona [RPM]

fprintf('Frequência elétrica nominal: f_e = %.1f Hz\n', f_e);
fprintf('Velocidade síncrona: w_sync = %.2f rad/s (%.0f RPM)\n', w_sync, N_sync);

% Corrente nominal aproximada (considerando fator de potência típico)
PF = 0.85; % Fator de potência típico
V_ll = V_dc / sqrt(3); % Tensão linha-linha aproximada
I_nom = P_nom / (sqrt(3) * V_ll * PF);
fprintf('Corrente nominal aproximada: I_nom = %.2f A\n', I_nom);

% Constante de tempo mecânica
T_m = J / B;
if B == 0
    fprintf('Constante de tempo mecânica: T_m = ∞ (B = 0)\n');
else
    fprintf('Constante de tempo mecânica: T_m = %.3f s\n', T_m);
end

T_ref_example = 10; % N·m
lambda_ref_example = 0.8; % Wb
[I_q_example, slip_example] = calculate_current_refs(T_ref_example, ...
    lambda_ref_example, P, Lm, L_rr);
fprintf('\nExemplo de cálculo de referências:\n');
fprintf('Para T_ref = %.1f N·m e λ_ref = %.1f Wb:\n', T_ref_example, lambda_ref_example);
fprintf('  I_q_ref = %.3f A\n', I_q_example);
fprintf('  slip = %.4f\n', slip_example);

%% FUNÇÃO AUXILIAR: CÁLCULO DE CORRENTES E TENSÕES
function [I_q_ref, slip] = calculate_current_refs(T_ref, lambda_ref, P, Lm, L_rr)
    % Calcula as correntes de referência para controle vetorial
    % I_d para fluxo
    I_d_ref = lambda_ref / Lm;
    
    % I_q para torque
    I_q_ref = T_ref / ((3/2) * (P/2) * (Lm/L_rr) * lambda_ref);
    
    % Escorregamento para controle indireto
    R_r = 0.816; % Valor fixo do problema
    slip = (R_r / L_rr) * (I_q_ref / I_d_ref);
end
