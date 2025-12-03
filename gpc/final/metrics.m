%% ═══════════════════════════════════════════════════════════════════════
%  SIMULAÇÃO COMPLETA: GPC dq → PWM abc → CORRENTES TRIFÁSICAS
%  Versão Corrigida e Otimizada
%% ═══════════════════════════════════════════════════════════════════════

clear; clc; close all;


%% PARÂMETROS DO SISTEMA
%*******************************************************
% Temporização
f = 60;                 % Frequência fundamental [Hz]
f_pwm = 10000;          % Frequência PWM [Hz]  
f_control = 1000;       % Frequência controle GPC [Hz]
Ts_pwm = 1/f_pwm;       % 100 µs
Ts_control = 1/f_control; % 1 ms
Tsim = 0.04;            % Tempo de simulação [s] - 40 ms

% Conversor
Vdc = 60;               % Tensão DC-link [V]
omega = 2*pi*f;         % Frequência angular [rad/s]

% Carga RL trifásica
R_load = 10;            % Resistência [Ω]
L_load = 0.010;         % Indutância [H] - 10 mH

fprintf('═══ CONFIGURAÇÃO ═══\n');
fprintf('PWM:       %.1f kHz (Ts = %.0f µs)\n', f_pwm/1e3, Ts_pwm*1e6);
fprintf('Controle:  %.1f kHz (Ts = %.0f µs)\n', f_control/1e3, Ts_control*1e6);
fprintf('Simulação: %.0f ms\n', Tsim*1e3);
fprintf('Carga:     R = %.1f Ω, L = %.1f mH\n', R_load, L_load*1e3);
fprintf('Vdc:       %.1f V\n\n', Vdc);

%% MATRIZES DE TRANSFORMAÇÃO
%*******************************************************
% Clarke inversa (αβ → abc)
K_clarke_inv = [1,        0;
                -1/2,     sqrt(3)/2;
                -1/2,    -sqrt(3)/2];

% Park (abc → dq) - para medição de correntes
% Usa transformação invariante em amplitude
K_park = @(theta) (2/3) * [cos(theta),   cos(theta - 2*pi/3),   cos(theta - 4*pi/3);
                           -sin(theta), -sin(theta - 2*pi/3), -sin(theta - 4*pi/3)];

%% VETORES DE TEMPO E INICIALIZAÇÃO
%*******************************************************
t_pwm = 0:Ts_pwm:Tsim;               % Vetor de tempo PWM
N_pwm = length(t_pwm);
t_control = 0:Ts_control:(Tsim-Ts_control);  % Vetor de tempo controle
N_control = length(t_control);

fprintf('Amostras PWM:      %d\n', N_pwm);
fprintf('Amostras controle: %d\n', N_control);

% ─────────────────────────────────────────────────────────────────────
% SINAIS DE CONTROLE dq (Saídas do GPC simuladas)
% ─────────────────────────────────────────────────────────────────────
% Aqui você substituiria pelo GPC real
% Por enquanto, sinais variantes no tempo para teste

Vd_gpc = 10 + 5*sin(2*pi*2*t_control);   % Vd variando [5, 15] V
Vq_gpc = 2 + 2*sin(2*pi*3*t_control);    % Vq variando [0, 4] V

fprintf('\nSinais GPC gerados: Vd ∈ [%.1f, %.1f] V, Vq ∈ [%.1f, %.1f] V\n', ...
        min(Vd_gpc), max(Vd_gpc), min(Vq_gpc), max(Vq_gpc));

% Arrays de saída
Va_ref = zeros(1, N_pwm);
Vb_ref = zeros(1, N_pwm);
Vc_ref = zeros(1, N_pwm);
ma = zeros(1, N_pwm);
mb = zeros(1, N_pwm);
mc = zeros(1, N_pwm);
PWM_a = zeros(1, N_pwm);
PWM_b = zeros(1, N_pwm);
PWM_c = zeros(1, N_pwm);

%% GERAÇÃO DA PORTADORA TRIANGULAR
%*******************************************************
% Gerar portadora diretamente em [0, 1] para evitar conversões
carrier = 0.5 + 0.5*sawtooth(2*pi*f_pwm*t_pwm, 0.5);  % Triangular [0,1]

fprintf('\n═══ INICIANDO CONVERSÃO dq → abc + PWM ═══\n');

%% LOOP PRINCIPAL: dq → αβ → abc → PWM
%*******************************************************
ratio = f_pwm / f_control;  % 10 amostras PWM por amostra GPC

for k_ctrl = 1:N_control
    
    % ─────────────────────────────────────────────────────────────────
    % FASE 1: Obter tensões de controle dq do GPC
    % ─────────────────────────────────────────────────────────────────
    Vd = Vd_gpc(k_ctrl);
    Vq = Vq_gpc(k_ctrl);
    
    % ─────────────────────────────────────────────────────────────────
    % FASE 2: Determinar índices PWM correspondentes
    % ─────────────────────────────────────────────────────────────────
    t_start = (k_ctrl - 1) * Ts_control;
    t_end = min(k_ctrl * Ts_control, Tsim);
    idx_pwm = find(t_pwm >= t_start & t_pwm < t_end);
    
    % ─────────────────────────────────────────────────────────────────
    % FASE 3: Para cada amostra PWM neste intervalo de controle
    % ─────────────────────────────────────────────────────────────────
    for jj = 1:length(idx_pwm)
        
        idx = idx_pwm(jj);
        t_current = t_pwm(idx);
        theta = omega * t_current;  % Ângulo elétrico instantâneo
        
        % ═════════════════════════════════════════════════════════════
        % TRANSFORMADA DE PARK INVERSA (dq → αβ)
        % ═════════════════════════════════════════════════════════════
        V_alpha = Vd * cos(theta) - Vq * sin(theta);
        V_beta  = Vd * sin(theta) + Vq * cos(theta);
        
        % ═════════════════════════════════════════════════════════════
        % TRANSFORMADA DE CLARKE INVERSA (αβ → abc)
        % ═════════════════════════════════════════════════════════════
        V_abc = K_clarke_inv * [V_alpha; V_beta];
        Va = V_abc(1);
        Vb = V_abc(2);
        Vc = V_abc(3);
        
        % Armazenar tensões de referência
        Va_ref(idx) = Va;
        Vb_ref(idx) = Vb;
        Vc_ref(idx) = Vc;
        
        % ═════════════════════════════════════════════════════════════
        % CÁLCULO DOS ÍNDICES DE MODULAÇÃO
        % ═════════════════════════════════════════════════════════════
        % CORREÇÃO CRÍTICA: Usar fórmula correta para SPWM bipolar
        % 
        % Para PWM bipolar com portadora unipolar [0,1]:
        % m = V_ref / V_dc + 0.5
        %
        % Onde V_ref ∈ [-V_dc, +V_dc] → m ∈ [0, 1]
        
        ma_val = Va / Vdc + 0.5;
        mb_val = Vb / Vdc + 0.5;
        mc_val = Vc / Vdc + 0.5;
        
        % Saturação rigorosa [0, 1]
        ma_val = max(min(ma_val, 1), 0);
        mb_val = max(min(mb_val, 1), 0);
        mc_val = max(min(mc_val, 1), 0);
        
        ma(idx) = ma_val;
        mb(idx) = mb_val;
        mc(idx) = mc_val;
        
        % ═════════════════════════════════════════════════════════════
        % COMPARAÇÃO COM PORTADORA → GERAÇÃO PWM
        % ═════════════════════════════════════════════════════════════
        carrier_val = carrier(idx);  % Já em [0, 1]
        
        PWM_a(idx) = ma_val > carrier_val;
        PWM_b(idx) = mb_val > carrier_val;
        PWM_c(idx) = mc_val > carrier_val;
    end
    
    % Progresso
    if mod(k_ctrl, floor(N_control/10)) == 0
        fprintf('  %.0f%% completo\n', (k_ctrl/N_control)*100);
    end
end

fprintf('✅ Conversão e PWM concluídos!\n');

%% VERIFICAÇÃO DOS ÍNDICES DE MODULAÇÃO
%*******************************************************
fprintf('\n═══ VERIFICAÇÃO DE ÍNDICES ═══\n');
fprintf('ma: [%.4f, %.4f]\n', min(ma), max(ma));
fprintf('mb: [%.4f, %.4f]\n', min(mb), max(mb));
fprintf('mc: [%.4f, %.4f]\n', min(mc), max(mc));

if all(ma >= 0 & ma <= 1) && all(mb >= 0 & mb <= 1) && all(mc >= 0 & mc <= 1)
    fprintf('✅ TODOS os índices em [0, 1]!\n');
else
    fprintf('❌ ERRO: Índices fora do range!\n');
    n_violations = sum(ma < 0 | ma > 1) + sum(mb < 0 | mb > 1) + sum(mc < 0 | mc > 1);
    fprintf('   Violações: %d de %d amostras\n', n_violations, 3*N_pwm);
end

%% SIMULAÇÃO DINÂMICA DAS CORRENTES (MODELO RL CORRETO)
%*******************************************************
fprintf('\n═══ SIMULANDO CORRENTES (Modelo RL) ═══\n');

% ─────────────────────────────────────────────────────────────────────
% TENSÕES APLICADAS PELA PONTE (Referência ao ponto médio DC)
% ─────────────────────────────────────────────────────────────────────
% Quando PWM = 1 → Chave superior ON → V_fase = +Vdc/2
% Quando PWM = 0 → Chave inferior ON → V_fase = -Vdc/2

Van = (PWM_a - 0.5) * Vdc;  % Tensão fase A-neutro
Vbn = (PWM_b - 0.5) * Vdc;
Vcn = (PWM_c - 0.5) * Vdc;

% ─────────────────────────────────────────────────────────────────────
% MODELO DIFERENCIAL: L*(dI/dt) = V - R*I
% ─────────────────────────────────────────────────────────────────────
% Resolver usando método de Euler para cada fase

Ia = zeros(1, N_pwm);
Ib = zeros(1, N_pwm);
Ic = zeros(1, N_pwm);

for k = 2:N_pwm
    % Equação: I[k+1] = I[k] + (Ts/L)*(V[k] - R*I[k])
    Ia(k) = Ia(k-1) + (Ts_pwm/L_load) * (Van(k-1) - R_load*Ia(k-1));
    Ib(k) = Ib(k-1) + (Ts_pwm/L_load) * (Vbn(k-1) - R_load*Ib(k-1));
    Ic(k) = Ic(k-1) + (Ts_pwm/L_load) * (Vcn(k-1) - R_load*Ic(k-1));
end

fprintf('✅ Correntes simuladas!\n');
fprintf('   Ia: [%.4f, %.4f] A\n', min(Ia), max(Ia));
fprintf('   Ib: [%.4f, %.4f] A\n', min(Ib), max(Ib));
fprintf('   Ic: [%.4f, %.4f] A\n', min(Ic), max(Ic));

% ─────────────────────────────────────────────────────────────────────
% FILTRAR CORRENTES PARA ANÁLISE (Remover ripple PWM)
% ─────────────────────────────────────────────────────────────────────
fc_filter = 500;  % Frequência de corte [Hz]
[b_filt, a_filt] = butter(4, 2*pi*fc_filter / (f_pwm*pi), 'low');

Ia_filt = filtfilt(b_filt, a_filt, Ia);  % filtfilt = zero phase
Ib_filt = filtfilt(b_filt, a_filt, Ib);
Ic_filt = filtfilt(b_filt, a_filt, Ic);

fprintf('   Ripple médio: %.4f A\n', mean(abs(Ia - Ia_filt)));

%% TRANSFORMADA PARA dq (MEDIÇÃO)
%*******************************************************
fprintf('\n═══ TRANSFORMANDO CORRENTES → dq ═══\n');

Id_meas = zeros(1, N_control);
Iq_meas = zeros(1, N_control);

for k_ctrl = 1:N_control
    t_sample = t_control(k_ctrl);
    
    % Encontrar amostra PWM mais próxima
    [~, idx_nearest] = min(abs(t_pwm - t_sample));
    
    theta_k = omega * t_sample;
    
    % Transformada de Park das correntes FILTRADAS
    Idq = K_park(theta_k) * [Ia_filt(idx_nearest); 
                              Ib_filt(idx_nearest); 
                              Ic_filt(idx_nearest)];
    
    Id_meas(k_ctrl) = Idq(1);
    Iq_meas(k_ctrl) = Idq(2);
end

fprintf('✅ Correntes dq calculadas!\n');
fprintf('   Id: [%.4f, %.4f] A\n', min(Id_meas), max(Id_meas));
fprintf('   Iq: [%.4f, %.4f] A\n', min(Iq_meas), max(Iq_meas));

%% REFERÊNCIAS TEÓRICAS (PARA COMPARAÇÃO)
%*******************************************************
% Em regime permanente senoidal: V = (R + jωL)*I
% Portanto: I = V / (R + jωL)

Z_magnitude = sqrt(R_load^2 + (omega*L_load)^2);
Z_phase = atan2(omega*L_load, R_load);

% Correntes esperadas (aproximação regime permanente)
Id_ref_theory = Vd_gpc / Z_magnitude;
Iq_ref_theory = Vq_gpc / Z_magnitude;

%% MÉTRICAS DE DESEMPENHO
%*******************************************************
fprintf('\n═══ MÉTRICAS DE DESEMPENHO ═══\n');

% THD (Total Harmonic Distortion) da corrente fase A
N_fft = 2^nextpow2(N_pwm);
Ia_fft = fft(Ia, N_fft);
f_axis = (0:N_fft-1) * f_pwm / N_fft;

% Fundamental (60 Hz)
[~, idx_fund] = min(abs(f_axis - f));
I_fundamental = abs(Ia_fft(idx_fund));

% Harmônicos
I_harmonics = Ia_fft;
I_harmonics(max(1, idx_fund-2):min(N_fft, idx_fund+2)) = 0;  % Zerar fundamental
I_harm_rms = sqrt(sum(abs(I_harmonics).^2) / N_fft);

THD = (I_harm_rms / I_fundamental) * 100;

fprintf('THD da corrente:    %.2f%%\n', THD);
fprintf('Corrente RMS (Ia):  %.4f A\n', rms(Ia));
fprintf('Ripple pico-pico:   %.4f A\n', max(Ia) - min(Ia));

% Erro de rastreamento em dq
erro_Id = rms(Id_meas - Id_ref_theory);
erro_Iq = rms(Iq_meas - Iq_ref_theory);

fprintf('Erro RMS Id:        %.4f A\n', erro_Id);
fprintf('Erro RMS Iq:        %.4f A\n', erro_Iq);

%% PLOTAGEM COMPLETA
%*******************************************************
fprintf('\n═══ GERANDO GRÁFICOS ═══\n');

figure('Position', [50, 50, 1600, 1200], 'Name', 'Simulação Completa GPC+PWM');

% ─────────────────────────────────────────────────────────────────────
% 1. SINAIS DE CONTROLE dq
% ─────────────────────────────────────────────────────────────────────
subplot(4,4,1);
plot(t_control*1e3, Vd_gpc, 'r-', 'LineWidth', 2); hold on;
plot(t_control*1e3, Vq_gpc, 'b-', 'LineWidth', 2);
title('Sinais de Controle GPC');
xlabel('Tempo [ms]'); ylabel('Tensão [V]');
legend('V_d', 'V_q', 'Location', 'best');
grid on;

% ─────────────────────────────────────────────────────────────────────
% 2. TENSÕES DE REFERÊNCIA abc
% ─────────────────────────────────────────────────────────────────────
subplot(4,4,2);
plot(t_pwm*1e3, Va_ref, 'r-', 'LineWidth', 0.5); hold on;
plot(t_pwm*1e3, Vb_ref, 'g-', 'LineWidth', 0.5);
plot(t_pwm*1e3, Vc_ref, 'b-', 'LineWidth', 0.5);
title('Tensões de Referência ABC');
xlabel('Tempo [ms]'); ylabel('Tensão [V]');
legend('V_a', 'V_b', 'V_c', 'Location', 'best');
grid on;

% ─────────────────────────────────────────────────────────────────────
% 3. ÍNDICES DE MODULAÇÃO
% ─────────────────────────────────────────────────────────────────────
subplot(4,4,3);
plot(t_pwm*1e3, ma, 'r-', 'LineWidth', 0.5); hold on;
plot(t_pwm*1e3, mb, 'g-', 'LineWidth', 0.5);
plot(t_pwm*1e3, mc, 'b-', 'LineWidth', 0.5);
plot([0 Tsim*1e3], [0 0], 'k--', 'LineWidth', 0.5);
plot([0 Tsim*1e3], [1 1], 'k--', 'LineWidth', 0.5);
title('Índices de Modulação');
xlabel('Tempo [ms]'); ylabel('Índice');
legend('m_a', 'm_b', 'm_c', 'Location', 'best');
grid on;
ylim([-0.1, 1.1]);

% ─────────────────────────────────────────────────────────────────────
% 4. HISTOGRAMA DOS ÍNDICES
% ─────────────────────────────────────────────────────────────────────
subplot(4,4,4);
histogram(ma, 30, 'FaceColor', 'r', 'FaceAlpha', 0.5); hold on;
histogram(mb, 30, 'FaceColor', 'g', 'FaceAlpha', 0.5);
histogram(mc, 30, 'FaceColor', 'b', 'FaceAlpha', 0.5);
xline(0, 'k--', 'LineWidth', 2);
xline(1, 'k--', 'LineWidth', 2);
title('Distribuição dos Índices');
xlabel('Índice'); ylabel('Frequência');
legend('m_a', 'm_b', 'm_c');
grid on;
xlim([-0.1, 1.1]);

% ─────────────────────────────────────────────────────────────────────
% 5. PWM vs PORTADORA (Zoom 2 ms)
% ─────────────────────────────────────────────────────────────────────
subplot(4,4,5);
t_zoom = t_pwm(t_pwm < 2e-3);
idx_zoom = 1:length(t_zoom);
plot(t_zoom*1e3, ma(idx_zoom), 'r-', 'LineWidth', 2); hold on;
plot(t_zoom*1e3, carrier(idx_zoom), 'k--', 'LineWidth', 1);
stairs(t_zoom*1e3, PWM_a(idx_zoom), 'b-', 'LineWidth', 1);
title('PWM Fase A (Zoom 2ms)');
xlabel('Tempo [ms]'); ylabel('Amplitude');
legend('Modulante', 'Portadora', 'PWM', 'Location', 'best');
grid on;
ylim([-0.1, 1.1]);

% ─────────────────────────────────────────────────────────────────────
% 6. SINAIS PWM (Zoom 5 ms)
% ─────────────────────────────────────────────────────────────────────
subplot(4,4,6);
t_zoom = t_pwm(t_pwm < 5e-3);
idx_zoom = 1:length(t_zoom);
stairs(t_zoom*1e3, PWM_a(idx_zoom), 'r-', 'LineWidth', 0.5); hold on;
stairs(t_zoom*1e3, PWM_b(idx_zoom) + 2, 'g-', 'LineWidth', 0.5);
stairs(t_zoom*1e3, PWM_c(idx_zoom) + 4, 'b-', 'LineWidth', 0.5);
title('Sinais PWM ABC (Zoom 5ms)');
xlabel('Tempo [ms]'); ylabel('Estado + Offset');
legend('PWM_A', 'PWM_B', 'PWM_C');
grid on;
ylim([-0.5, 5]);

% ─────────────────────────────────────────────────────────────────────
% 7. CORRENTES INSTANTÂNEAS (Zoom 10 ms)
% ─────────────────────────────────────────────────────────────────────
subplot(4,4,7);
t_zoom = t_pwm(t_pwm < 10e-3);
idx_zoom = 1:length(t_zoom);
plot(t_zoom*1e3, Ia(idx_zoom), 'r-', 'LineWidth', 0.5); hold on;
plot(t_zoom*1e3, Ib(idx_zoom), 'g-', 'LineWidth', 0.5);
plot(t_zoom*1e3, Ic(idx_zoom), 'b-', 'LineWidth', 0.5);
title('Correntes ABC (Zoom 10ms)');
xlabel('Tempo [ms]'); ylabel('Corrente [A]');
legend('I_a', 'I_b', 'I_c', 'Location', 'best');
grid on;

% ─────────────────────────────────────────────────────────────────────
% 8. CORRENTES FILTRADAS (Todo o tempo)
% ─────────────────────────────────────────────────────────────────────
subplot(4,4,8);
plot(t_pwm*1e3, Ia_filt, 'r-', 'LineWidth', 1.5); hold on;
plot(t_pwm*1e3, Ib_filt, 'g-', 'LineWidth', 1.5);
plot(t_pwm*1e3, Ic_filt, 'b-', 'LineWidth', 1.5);
title('Correntes Filtradas (Fundamental)');
xlabel('Tempo [ms]'); ylabel('Corrente [A]');
legend('I_a', 'I_b', 'I_c', 'Location', 'best');
grid on;

% ─────────────────────────────────────────────────────────────────────
% 9. CORRENTE vs TENSÃO (Fase A)
% ─────────────────────────────────────────────────────────────────────
subplot(4,4,9);
t_zoom = t_pwm(t_pwm < 20e-3);
idx_zoom = 1:length(t_zoom);
yyaxis left;
plot(t_zoom*1e3, Ia_filt(idx_zoom), 'r-', 'LineWidth', 2);
ylabel('Corrente I_a [A]');
yyaxis right;
plot(t_zoom*1e3, Van(idx_zoom), 'b-', 'LineWidth', 0.5);
ylabel('Tensão V_{an} [V]');
title('Tensão vs Corrente (Fase A)');
xlabel('Tempo [ms]');
grid on;

% ─────────────────────────────────────────────────────────────────────
% 10. ESPECTRO DA CORRENTE FASE A
% ─────────────────────────────────────────────────────────────────────
subplot(4,4,10);
plot(f_axis(1:N_fft/2), 20*log10(abs(Ia_fft(1:N_fft/2))), 'b-', 'LineWidth', 1);
hold on;
xline(f, 'r--', 'LineWidth', 2, 'Label', 'Fundamental');
xline(f_pwm, 'k--', 'LineWidth', 1, 'Label', 'f_{PWM}');
title(sprintf('Espectro I_a (THD = %.2f%%)', THD));
xlabel('Frequência [Hz]'); ylabel('Magnitude [dB]');
xlim([0, 3000]);
grid on;

% ─────────────────────────────────────────────────────────────────────
% 11. CORRENTES dq MEDIDAS
% ─────────────────────────────────────────────────────────────────────
subplot(4,4,11);
plot(t_control*1e3, Id_meas, 'r-', 'LineWidth', 2); hold on;
plot(t_control*1e3, Iq_meas, 'b-', 'LineWidth', 2);
title('Correntes Medidas (dq)');
xlabel('Tempo [ms]'); ylabel('Corrente [A]');
legend('I_d', 'I_q', 'Location', 'best');
grid on;

% ─────────────
