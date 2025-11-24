clear; clc; close all;

%% Parâmetros do PMSM  - Scooter Peugeot
P_nom = 11e3;       % Potencia nominal W
Rs = 0.5;           % Resistência do estator [Ω]
Ld = 20.1e-3;       % Indutância do eixo d [H]
Lq = 40.9e-3;       % Indutância do eixo q [H] 
lambda_m = 0.5126;  % Fluxo magnetico concatenado do ímã permanente [V/rad/s]
lambda_m_rpm = 278.92;  % Fluxo em rpm V.s/krpm;
psi_PM_rpm = lambda_m_rpm ; 

P = 6;              % Número de polos
lambda_pm = 60 / (sqrt(3)*pi*P*1e3) * psi_PM_rpm ;

J = 0.03877;        % Inércia do rotor [kg.m²]
B = 0.0194;         % Coeficiente de atrito [N.m.s]

% Parâmetros do sistema de controle
V_battery = 96;     % Tensão da bateria [V]
V_in = 5;           % Tensão de entrada do modulador [V]
K_PWM = 1/V_in;     % Ganho do modulador PWM
K_r = V_battery;    % Ganho do inversor (tensão CC)


%% Função de transferência da planta de corrente (eixo q)
K_iq = 1;           % Ganho do sensor de corrente

% G_iq(s) = 1 / (Lq*s + Rs)
num_G_iq = 1;
den_G_iq = [Lq, Rs];
G_iq = tf(num_G_iq, den_G_iq);

%% Função de transferência em laço aberto não compensada
% FTLA_nc(s) = K_PWM * K_r * G_iq(s) * K_iq
FTLA_nc_iq = K_PWM * K_r * G_iq * K_iq;

disp('Função de transferência em laço aberto não compensada:');
FTLA_nc_iq


%% Projeto do compensador de corrente C_iq(s)
% Usando um controlador PI para melhorar o desempenho
% C_iq(s) = Kp + Ki/s = (Kp*s + Ki)/s

% Dados
fs = 20000;                % Frequência de chaveamento (Hz)
fciq = fs/10;              % Frequência de corte (Hz)
wciq = 2*pi*fciq;          % Velocidade angular de corte (rad/s)
Mphi = pi/3;               % Margem de fase (rad)

tau_r = 1/2 * 1/fs;

% Especificações de projeto:
BW_current = fciq;  % Largura de banda desejada [Hz]
PM_desired = rad2deg(Mphi);    % Margem de fase desejada [graus]

% % Cálculo dos ganhos do PI
% wc = 2*pi*BW_current;  % Frequência de cruzamento [rad/s]
% wc_current_iq=wc;
% FTLA_nc_iq_at_wc = sqrt(-1); % Valor de ejemplo

% [mag, phase] = bode(FTLA_nc, angle(FTLA_nc_iq_at_wc)*wciq);
% phase_rad = phase * pi/180; % fase em rad/s
% wziq = wc_current_iq / tan(Mphi - pi/2 - phase_rad)
% kc_iq = wciq / (sqrt(wciq^2 + wziq^2) * angle(FTLA_nc_iq_at_wc)*mag);
% k_c_iq = wciq / (sqrt(wciq^2 + wziq^2) * angle(FTLA_nc_iq_at_wc));

%%
% Obter fase na frequência de cruzamento
[mag, phase] = bode(FTLA_nc_iq, wciq);
mag_wc = mag;
mag_at_wc = abs(freqresp(FTLA_nc_iq, wciq)); % outro jeito de calcular magnitude
phase_rad = phase * pi/180; % fase em rad/s

wziq = wciq / tan(Mphi - pi/2 - phase_rad);
disp('wziq rad/s:');
wziq

tau_wziq = 1 / wziq;

FTLA_nc_iq_at_wc = sqrt(-1); 

% Cálculo do ganho kc
k_c_iq = wciq / (sqrt(wciq^2 + wziq^2) * mag_wc);
k_c_iq
kc_iq = wciq / (sqrt(wciq^2 + wziq^2) * angle(FTLA_nc_iq_at_wc)*mag)
kc_iq
% Ganho necessário na frequência de cruzamento
Kp = 1/mag_at_wc;

% Zero do PI em wc/5 para boa margem de fase
Ki = Kp * (wciq/5);

% Compensador PI
% i_qs / v_qs 
C_iq = tf([Kp, Ki], [1, 0]);

disp('Compensador de corrente C_iq(s):');
C_iq
C_iq_zpk = zpk(C_iq);
disp('Forma com zpk:');
C_iq_zpk

% Função de transferência em laço aberto compensada
% FTLA_c(s) = C_iq(s) * K_PWM * K_r * G_iq(s) * K_iq
FTLA_c_iq = C_iq * K_PWM * K_r * G_iq * K_iq;

disp('Função de transferência em laço aberto compensada:');
FTLA_c_iq

%% Análise no domínio da frequência
figure('Position', [100, 100, 1200, 800]);

% Diagrama de Bode
subplot(2,2,1);
bode(FTLA_nc_iq, 'b', FTLA_c_iq, 'r', {1, 1e5});
grid on;
legend('Não Compensada', 'Compensada', 'Location', 'best');
title('Diagrama de Bode - FTLA');

% Margens de estabilidade
subplot(2,2,2);
margin(FTLA_nc_iq);
grid on;
title('Margens de Estabilidade - Não Compensada');

subplot(2,2,3);
margin(FTLA_c_iq);
grid on;
title('Margens de Estabilidade - Compensada');

% Resposta ao degrau em malha fechada
subplot(2,2,4);
FTMF_nc_iq = feedback(FTLA_nc_iq, 1);
FTMF_c_iq = feedback(FTLA_c_iq, 1);
step(FTMF_nc_iq, 'b', FTMF_c_iq, 'r', 0.1);
grid on;
legend('Não Compensada', 'Compensada', 'Location', 'best');
title('Resposta ao Degrau em Malha Fechada');
xlabel('Tempo [s]');
ylabel('Amplitude');

%%******************************************************************************
%%
%%******************************************************************************
%% Função de transferência da planta de corrente (eixo d)
K_id = 1;           % Ganho do sensor de corrente

% G_iq(s) = 1 / (Lq*s + Rs)
num_G_id = 1;
den_G_id = [Ld, Rs];
G_id = tf(num_G_id, den_G_id);
FTLA_nc_id = K_PWM * K_r * G_id * K_id;

fcid = fs/10;              % Frequência de corte (Hz)
wcid = 2*pi*fcid;          % Velocidade angular de corte (rad/s)
Mphi = pi/3;               % Margem de fase (rad)
BW_current = fcid;  % Largura de banda desejada [Hz]
PM_desired = rad2deg(Mphi);    % Margem de fase desejada [graus]


% Obter fase na frequência de cruzamento
[mag, phase] = bode(FTLA_nc_id, wcid);
mag_wc = mag;
mag_at_wc = abs(freqresp(FTLA_nc_id, wcid)); % outro jeito de calcular magnitude
phase_rad = phase * pi/180; % fase em rad/s

wzid = wcid / tan(Mphi - pi/2 - phase_rad);
disp('wzid rad/s:');
wzid
tau_wzid = 1 / wzid;

FTLA_nc_id_at_wc = sqrt(-1); 
% Cálculo do ganho kc
k_c_id = wcid / (sqrt(wcid^2 + wzid^2) * mag_wc);
k_c_id
kc_id = wcid / (sqrt(wcid^2 + wzid^2) * angle(FTLA_nc_id_at_wc)*mag)
kc_id

% Ganho necessário na frequência de cruzamento
Kp = 1/mag_at_wc;
% Zero do PI em wc/5 para boa margem de fase
Ki = Kp * (wcid/5);

% Compensador PI
% i_qs / v_qs 
C_id = tf([Kp, Ki], [1, 0]);

disp('Compensador de corrente C_id(s):');
C_id

% Função de transferência em laço aberto compensada
% FTLA_c(s) = C_id(s) * K_PWM * K_r * G_iq(s) * K_id
FTLA_c_id = C_id * K_PWM * K_r * G_id * K_id;

disp('Função de transferência em laço aberto compensada:');
FTLA_c_id
FTLA_c_id_zpk = zpk(FTLA_c_id);
disp('Forma com zpk:');
FTLA_c_id_zpk



%%******************************************************************************
%%
% Planta de velocidade
%%******************************************************************************
%% Planta de velocidade G_ni(s) = nrpm / iqs
% G_ni = [30/pi * 3/2 * P/2 * lambda_m] * 1/(J*s)

K_vel = (30/pi) * (3/2) * (P/2) * lambda_m;
num_G_ni = K_vel;
den_G_ni = [J, 0];  % 1/(J*s)
G_ni  = tf(num_G_ni, den_G_ni);

disp('Planta de velocidade G_ni(s):');
G_ni
G_ni_zpk = zpk(G_ni);
disp('Forma com zpk:');
G_ni_zpk

% Ganho do sensor de velocidade K_ni
% Assumindo sensor que converte rpm para tensão (ex: 10V = 3000 rpm)
K_ni = 1;  % [V/rpm] - ajustável
f_cn = fciq/10;     % Frequência de corte da malha de velocidade 50x menor que corrente f. [Hz]
w_cn = 2*pi*f_cn;   % Frequência de corte [rad/s]
M_Phi = pi/3;       % Margem de fase desejada [rad]
omega_cn = M_Phi;

fprintf('\n ESPECIFICAÇÕES DO CONTROLE DE VELOCIDADE \n');
fprintf('Frequência de corte: f_cn = %.1f Hz\n', f_cn);
fprintf('Frequência angular de corte: w_cn = %.2f rad/s\n', w_cn);
fprintf('Margem de fase desejada: M_Phi = %.2f°, rad/s\n', M_Phi*180/pi, M_Phi);

FTLA_nc_N = (1/K_iq) * G_ni_zpk * K_ni;

% Cálculo do compensador de velocidade
% Ângulo da FTLA_nc na frequência de corte
[mag_wcn, phase_wcn] = bode(FTLA_nc_N, w_cn);
phase_wcn_rad = phase_wcn * pi/180;  % Converter para radianos

fprintf('\nÂngulo da FTLA_nc em w_cn: %.2f°\n', phase_wcn);

% Cálculo do zero do compensador
w_zn = w_cn / tan(M_Phi - pi/2 - phase_wcn_rad);

% Cálculo do ganho do compensador
K_cn = (w_cn / sqrt(w_cn^2 + w_zn^2)) * (1 / abs(mag_wcn));

% Compensador PI
% C_n(s) = K_cn * (s + w_zn) / s
C_n = tf(K_cn * [1, w_zn], [1, 0]);

disp('Compensador de velocidade C_n(s) e com formato zpk:');
C_n
zpk(C_n)

% Obter os coeficientes do numerador e denominador
[num, den] = tfdata(C_n, 'v');

Kp = num(1);  % Coeficiente de s no numerador
Ki = num(2);  % Termo constante no numerador

fprintf('\n PARÂMETROS DO CONTROLADOR PI Velocidade\n');
fprintf('Ganho Proporcional (Kp): %.4f\n', Kp);
fprintf('Ganho Integral (Ki): %.4f\n', Ki);

% Função de transferência em laço aberto compensada - Velocidade
FTLA_c_N = C_n * FTLA_nc_N;

disp('Função de transferência em laço aberto compensada (velocidade):');
FTLA_c_N

%% Análise no domínio da frequência - Controle de Velocidade
figure('Position', [100, 100, 1200, 800]);

% Diagrama de Bode - Velocidade
subplot(2,2,1);
bode(FTLA_nc_N, 'b', FTLA_c_N, 'r', {1, 1000});
grid on;
legend('Não Compensada', 'Compensada', 'Location', 'best');
title('Diagrama de Bode - Controle de Velocidade');

% Margens de estabilidade - Velocidade
subplot(2,2,2);
margin(FTLA_nc_N);
grid on;
title('Margens - Não Compensada (Velocidade)');

subplot(2,2,3);
margin(FTLA_c_N);
grid on;
title('Margens - Compensada (Velocidade)');

% Resposta ao degrau em malha fechada - Velocidade
subplot(2,2,4);
FTMF_nc_N = feedback(FTLA_nc_N, 1);
FTMF_c_N = feedback(FTLA_c_N, 1);
step(FTMF_nc_N, 'b', FTMF_c_N, 'r', 2);
grid on;
legend('Não Compensada', 'Compensada', 'Location', 'best');
title('Resposta ao Degrau - Controle de Velocidade');
xlabel('Tempo [s]');
ylabel('Velocidade [rpm]');

% Análise detalhada das margens de estabilidade - Velocidade
[Gm_nc_N, Pm_nc_N, Wcg_nc_N, Wcp_nc_N] = margin(FTLA_nc_N);
[Gm_c_N, Pm_c_N, Wcg_c_N, Wcp_c_N] = margin(FTLA_c_N);

fprintf('\n=== ANÁLISE DE ESTABILIDADE - CONTROLE DE VELOCIDADE ===\n');
fprintf('Sistema NÃO COMPENSADO:\n');
fprintf('  Margem de Ganho: %.2f dB @ %.2f rad/s\n', 20*log10(Gm_nc_N), Wcg_nc_N);
fprintf('  Margem de Fase: %.2f° @ %.2f rad/s\n', Pm_nc_N, Wcp_nc_N);

fprintf('\nSistema COMPENSADO:\n');
fprintf('  Margem de Ganho: %.2f dB @ %.2f rad/s\n', 20*log10(Gm_c_N), Wcg_c_N);
fprintf('  Margem de Fase: %.2f° @ %.2f rad/s\n', Pm_c_N, Wcp_c_N);

