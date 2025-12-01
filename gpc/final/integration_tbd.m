%% CÓDIGO COMPLETO CORRIGIDO - CONVERSÃO + PWM + CORRENTES
clear; clc; close all;

% Parâmetros do sistema
f = 60;              % Frequência fundamental [Hz]
f_pwm = 10000;       % Frequência PWM [Hz]  
f_control = 1000;    % Frequência controle [Hz]
Ts_pwm = 1/f_pwm;
Ts_control = 1/f_control;
Tsim = 0.04;         % Tempo de simulação [s]

Vdc = 60;            % Tensão DC-link [V]
omega = 2*pi*f;      % Frequência angular

% Matriz de Clarke inversa
K_inv = [1, 0; -1/2, sqrt(3)/2; -1/2, -sqrt(3)/2];

% Vetores de tempo
t = 0:Ts_pwm:Tsim;
N = length(t);

% Inicializar arrays
u_d_array = 10 + 3*sin(2*pi*2*t(1:floor(Tsim/Ts_control))); % Exemplo Vd
u_q_array = 2 + 1*sin(2*pi*3*t(1:floor(Tsim/Ts_control)));  % Exemplo Vq

va_ref = zeros(1,N); vb_ref = zeros(1,N); vc_ref = zeros(1,N);
ma = zeros(1,N); mb = zeros(1,N); mc = zeros(1,N);
ga = zeros(1,N); gb = zeros(1,N); gc = zeros(1,N);

% Portadora triangular
carrier = sawtooth(2*pi*f_pwm*t + pi/2, 0.5);

fprintf('=== INICIANDO SIMULAÇÃO ===\n');
fprintf('Tempo total: %.0f ms\n', Tsim*1000);
fprintf('Amostras PWM: %d\n', N);
fprintf('Amostras controle: %d\n', length(u_d_array));

%% CONVERSÃO + PWM (CÓDIGO CORRIGIDO)
ind = length(u_d_array);

for k = 1:ind
    u_d = u_d_array(k);
    u_q = u_q_array(k);
    t_start = (k-1)*Ts_control;
    t_end = min(k*Ts_control, Tsim);
    idxs = find(t >= t_start & t < t_end);
    
    for jj = 1:length(idxs)
        ti = t(idxs(jj));
        theta = omega * ti;
        
        % 1. Inverse Park (dq → αβ)
        v_alpha = cos(theta)*u_d - sin(theta)*u_q;
        v_beta  = sin(theta)*u_d + cos(theta)*u_q;
        
        % 2. Inverse Clarke (αβ → abc)
        vabc = K_inv * [v_alpha; v_beta];
        va = vabc(1); vb = vabc(2); vc = vabc(3);
        
        % 3. CORREÇÃO CRÍTICA: Normalização para modulação 0-1
        % Para SPWM com offset, precisamos garantir 0 ≤ m ≤ 1
        ma_i = (va / (Vdc/2) + 1) / 2;  % Normaliza para [0,1]
        mb_i = (vb / (Vdc/2) + 1) / 2;
        mc_i = (vc / (Vdc/2) + 1) / 2;
        
        % 4. Saturação rigorosa entre 0 e 1
        ma_i = max(min(ma_i, 1), 0);
        mb_i = max(min(mb_i, 1), 0); 
        mc_i = max(min(mc_i, 1), 0);
        
        % 5. Armazenar valores
        va_ref(idxs(jj)) = va;
        vb_ref(idxs(jj)) = vb; 
        vc_ref(idxs(jj)) = vc;
        ma(idxs(jj)) = ma_i;
        mb(idxs(jj)) = mb_i;
        mc(idxs(jj)) = mc_i;
        
        % 6. Geração PWM (comparação com portadora 0-1)
        carrier_val = (carrier(idxs(jj)) + 1) / 2;  % Portadora em [0,1]
        ga(idxs(jj)) = ma_i > carrier_val;
        gb(idxs(jj)) = mb_i > carrier_val;
        gc(idxs(jj)) = mc_i > carrier_val;
    end
    
    if mod(k, 10) == 0
        fprintf('Processado: %.1f%%\n', (k/ind)*100);
    end
end

fprintf('Conversão e PWM concluídos!\n');

%% SIMULAÇÃO DAS CORRENTES
fprintf('\n=== SIMULANDO CORRENTES ===\n');

% Parâmetros da carga
R_load = 10;     % Resistência [Ω]
L_load = 0.01;   % Indutância [H]

% Tensões de fase do conversor
V_an = (ga - 1/3) * Vdc;  % Tensão fase A-neutro
V_bn = (gb - 1/3) * Vdc;  % Tensão fase B-neutro  
V_cn = (gc - 1/3) * Vdc;  % Tensão fase C-neutro

% Simulação das correntes (modelo simplificado)
I_a = zeros(1,N);
I_b = zeros(1,N); 
I_c = zeros(1,N);

% Filtro para suavizar as correntes
[b, a] = butter(2, 2*pi*500/(f_pwm/2));

I_a_filtered = filter(b, a, V_an) / R_load;
I_b_filtered = filter(b, a, V_bn) / R_load;
I_c_filtered = filter(b, a, V_cn) / R_load;

fprintf('Correntes simuladas!\n');

%% VERIFICAÇÃO DOS ÍNDICES DE MODULAÇÃO
fprintf('\n=== VERIFICAÇÃO DOS ÍNDICES DE MODULAÇÃO ===\n');

fprintf('Índices de modulação:\n');
fprintf('ma: [%.3f, %.3f] ✓\n', min(ma), max(ma));
fprintf('mb: [%.3f, %.3f] ✓\n', min(mb), max(mb)); 
fprintf('mc: [%.3f, %.3f] ✓\n', min(mc), max(mc));

if all(ma >= 0 & ma <= 1) && all(mb >= 0 & mb <= 1) && all(mc >= 0 & mc <= 1)
    fprintf('✅ TODOS os índices entre 0 e 1!\n');
else
    fprintf('❌ ALGUM índice fora do range [0,1]!\n');
end

%% PLOTAGEM COMPLETA DOS RESULTADOS
fprintf('\n=== GERANDO GRÁFICOS ===\n');

figure('Position', [50, 50, 1400, 1000]);

% 1. Sinais de controle dq
subplot(4,3,1);
t_control_plot = 0:Ts_control:(length(u_d_array)-1)*Ts_control;
plot(t_control_plot, u_d_array, 'r-', 'LineWidth', 2); hold on;
plot(t_control_plot, u_q_array, 'b-', 'LineWidth', 2);
title('Sinais de Controle dq');
xlabel('Tempo [s]');
ylabel('Tensão [V]');
legend('V_d', 'V_q', 'Location', 'best');
grid on;

% 2. Tensões de referência abc
subplot(4,3,2);
plot(t, va_ref, 'r-', 'LineWidth', 1); hold on;
plot(t, vb_ref, 'g-', 'LineWidth', 1);
plot(t, vc_ref, 'b-', 'LineWidth', 1);
title('Tensões de Referência ABC');
xlabel('Tempo [s]');
ylabel('Tensão [V]');
legend('V_a', 'V_b', 'V_c', 'Location', 'best');
grid on;

% 3. Índices de modulação
subplot(4,3,3);
plot(t, ma, 'r-', 'LineWidth', 1); hold on;
plot(t, mb, 'g-', 'LineWidth', 1);
plot(t, mc, 'b-', 'LineWidth', 1);
title('Índices de Modulação ABC');
xlabel('Tempo [s]');
ylabel('Índice de Modulação');
legend('m_a', 'm_b', 'm_c', 'Location', 'best');
grid on;
ylim([-0.1, 1.1]);

% 4. Verificação dos índices (histograma)
subplot(4,3,4);
histogram(ma, 20, 'FaceColor', 'r', 'FaceAlpha', 0.7); hold on;
histogram(mb, 20, 'FaceColor', 'g', 'FaceAlpha', 0.7);
histogram(mc, 20, 'FaceColor', 'b', 'FaceAlpha', 0.7);
title('Distribuição dos Índices de Modulação');
xlabel('Índice de Modulação');
ylabel('Frequência');
legend('m_a', 'm_b', 'm_c', 'Location', 'best');
grid on;
xlim([-0.1, 1.1]);

% 5. PWM e portadora (zoom)
subplot(4,3,5);
t_zoom = t(t < 0.002);
idx_zoom = 1:length(t_zoom);
plot(t_zoom, ma(idx_zoom), 'r-', 'LineWidth', 2); hold on;
plot(t_zoom, (carrier(idx_zoom)+1)/2, 'k--', 'LineWidth', 1);
plot(t_zoom, ga(idx_zoom), 'b-', 'LineWidth', 1);
title('PWM - Modulante vs Portadora (Zoom)');
xlabel('Tempo [s]');
ylabel('Amplitude');
legend('m_a', 'Portadora', 'PWM_A', 'Location', 'best');
grid on;

% 6. Sinais PWM
subplot(4,3,6);
t_zoom = t(t < 0.005);
idx_zoom = 1:length(t_zoom);
plot(t_zoom, ga(idx_zoom), 'r-', 'LineWidth', 1); hold on;
plot(t_zoom, gb(idx_zoom), 'g-', 'LineWidth', 1);
plot(t_zoom, gc(idx_zoom), 'b-', 'LineWidth', 1);
title('Sinais PWM Gerados (Zoom)');
xlabel('Tempo [s]');
ylabel('Estado (0/1)');
legend('PWM_A', 'PWM_B', 'PWM_C', 'Location', 'best');
grid on;
ylim([-0.1, 1.1]);

% 7. Correntes de saída (zoom inicial)
subplot(4,3,7);
t_zoom = t(t < 0.01);
idx_zoom = 1:length(t_zoom);
plot(t_zoom, I_a_filtered(idx_zoom), 'r-', 'LineWidth', 1); hold on;
plot(t_zoom, I_b_filtered(idx_zoom), 'g-', 'LineWidth', 1);
plot(t_zoom, I_c_filtered(idx_zoom), 'b-', 'LineWidth', 1);
title('Correntes de Saída (Zoom Inicial)');
xlabel('Tempo [s]');
ylabel('Corrente [A]');
legend('I_a', 'I_b', 'I_c', 'Location', 'best');
grid on;

% 8. Correntes de saída (todo o tempo)
subplot(4,3,8);
plot(t, I_a_filtered, 'r-', 'LineWidth', 1); hold on;
plot(t, I_b_filtered, 'g-', 'LineWidth', 1);
plot(t, I_c_filtered, 'b-', 'LineWidth', 1);
title('Correntes de Saída (Tempo Total)');
xlabel('Tempo [s]');
ylabel('Corrente [A]');
legend('I_a', 'I_b', 'I_c', 'Location', 'best');
grid on;

% 9. Espectro das correntes
subplot(4,3,9);
N_fft = 2^nextpow2(length(I_a_filtered));
f_axis = (0:N_fft-1) * f_pwm / N_fft;
I_a_fft = fft(I_a_filtered, N_fft);
plot(f_axis(1:N_fft/2), abs(I_a_fft(1:N_fft/2)), 'r-', 'LineWidth', 1);
title('Espectro da Corrente Fase A');
xlabel('Frequência [Hz]');
ylabel('Magnitude');
xlim([0, 2000]);
grid on;

% 10. Transformada das correntes para dq
subplot(4,3,10);
Id_measured = zeros(1, ind);
Iq_measured = zeros(1, ind);

for k = 1:ind
    t_mid = (k-0.5)*Ts_control;
    idx_mid = find(t >= t_mid, 1);
    if ~isempty(idx_mid) && idx_mid <= length(I_a_filtered)
        theta_k = omega * t_mid;
        % Transformada Park das correntes medidas
        Id_measured(k) = (2/3) * (I_a_filtered(idx_mid)*cos(theta_k) + ...
                                 I_b_filtered(idx_mid)*cos(theta_k - 2*pi/3) + ...
                                 I_c_filtered(idx_mid)*cos(theta_k - 4*pi/3));
        Iq_measured(k) = (2/3) * (-I_a_filtered(idx_mid)*sin(theta_k) - ...
                                 I_b_filtered(idx_mid)*sin(theta_k - 2*pi/3) - ...
                                 I_c_filtered(idx_mid)*sin(theta_k - 4*pi/3));
    end
end

plot(t_control_plot, Id_measured, 'r-', 'LineWidth', 2); hold on;
plot(t_control_plot, Iq_measured, 'b-', 'LineWidth', 2);
title('Correntes Medidas nos Eixos dq');
xlabel('Tempo [s]');
ylabel('Corrente [A]');
legend('I_d medido', 'I_q medido', 'Location', 'best');
grid on;

% 11. Comparação referência vs medida
subplot(4,3,11);
% Criar referências de corrente baseadas nas tensões de controle
Id_ref = u_d_array / R_load;  % Simplificação
Iq_ref = u_q_array / R_load;

plot(t_control_plot, Id_ref, 'r--', 'LineWidth', 2); hold on;
plot(t_control_plot, Id_measured, 'r-', 'LineWidth', 1);
plot(t_control_plot, Iq_ref, 'b--', 'LineWidth', 2);
plot(t_control_plot, Iq_measured, 'b-', 'LineWidth', 1);
title('Comparação: Referência vs Medida');
xlabel('Tempo [s]');
ylabel('Corrente [A]');
legend('I_d ref', 'I_d med', 'I_q ref', 'I_q med', 'Location', 'best');
grid on;

% 12. Resumo
subplot(4,3,12);
text(0.1, 0.9, 'RESUMO DA SIMULAÇÃO', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.7, sprintf('Tempo: %.0f ms', Tsim*1000));
text(0.1, 0.6, sprintf('PWM: %.1f kHz', f_pwm/1000));
text(0.1, 0.5, sprintf('Controle: %.0f Hz', f_control));
text(0.1, 0.4, sprintf('m_a: [%.3f, %.3f]', min(ma), max(ma)));
text(0.1, 0.3, sprintf('m_b: [%.3f, %.3f]', min(mb), max(mb)));
text(0.1, 0.2, sprintf('m_c: [%.3f, %.3f]', min(mc), max(mc)));
text(0.1, 0.1, sprintf('I_max: %.2f A', max([I_a_filtered, I_b_filtered, I_c_filtered])));
axis off;

fprintf('\n✅ SIMULAÇÃO CONCLUÍDA COM SUCESSO!\n');
