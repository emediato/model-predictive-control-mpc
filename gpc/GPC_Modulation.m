%% Geração dos Sinais de Referência e Portadora
t = 0:Ts:Tsim;
N = length(t);

% Índice de modulação e frequência
m = 0.9;                % Índice de modulação (0 < m <= 1)
f_ref = f_fundamental;  % Frequência da referência

% Sinais de referência trifásicos
Vref_A = m * sin(2*pi*f_ref*t);
Vref_B = m * sin(2*pi*f_ref*t - 2*pi/3);
Vref_C = m * sin(2*pi*f_ref*t + 2*pi/3);

% Portadora triangular
fprintf('Gerando portadora triangular...\n');
carrier = sawtooth(2*pi*f_chaveamento*t, 0.5); % Triangular centrada em 0

% Portadoras deslocadas para modulação PWM
carrier_upper = carrier;           % Portadora superior
carrier_lower = -carrier;          % Portadora inferior

%% Técnicas de Modulação PWM
fprintf('Aplicando técnicas de modulação...\n');

% 1. PWM Senoidal Tradicional
fprintf('1. PWM Senoidal Tradicional\n');
PWM_A_SPWM = zeros(size(t));
PWM_B_SPWM = zeros(size(t));
PWM_C_SPWM = zeros(size(t));

for k = 1:length(t)
    if Vref_A(k) > carrier(k)
        PWM_A_SPWM(k) = 1;
    else
        PWM_A_SPWM(k) = 0;
    end
    
    if Vref_B(k) > carrier(k)
        PWM_B_SPWM(k) = 1;
    else
        PWM_B_SPWM(k) = 0;
    end
    
    if Vref_C(k) > carrier(k)
        PWM_C_SPWM(k) = 1;
    else
        PWM_C_SPWM(k) = 0;
    end
end

% 2. PWM Vetorial Space Vector PWM (SVPWM)
fprintf('2. Space Vector PWM (SVPWM)\n');
% Cálculo das tensões de referência no plano αβ
V_alpha = (2/3) * (Vref_A - 0.5*Vref_B - 0.5*Vref_C);
V_beta = (2/3) * (sqrt(3)/2*Vref_B - sqrt(3)/2*Vref_C);

% Normalização para o hexágono
Vref_mag = sqrt(V_alpha.^2 + V_beta.^2);
Vref_angle = atan2(V_beta, V_alpha);

% Implementação simplificada do SVPWM
[PWM_A_SVPWM, PWM_B_SVPWM, PWM_C_SVPWM] = svpwm_simplified(V_alpha, V_beta, carrier, Vdc);

% 3. PWM Terceira Harmônica Injetada (THIPWM)
fprintf('3. THIPWM - Terceira Harmônica Injetada\n');
third_harmonic = 0.25 * sin(3*2*pi*f_ref*t); % Injeção de 3ª harmônica
Vref_A_THI = Vref_A + third_harmonic;
Vref_B_THI = Vref_B + third_harmonic;
Vref_C_THI = Vref_C + third_harmonic;

% Normalização
max_ref = max([max(Vref_A_THI), max(Vref_B_THI), max(Vref_C_THI)]);
min_ref = min([min(Vref_A_THI), min(Vref_B_THI), min(Vref_C_THI)]);
scale_factor = 2/(max_ref - min_ref);
offset = (max_ref + min_ref)/2;

Vref_A_THI = (Vref_A_THI - offset) * scale_factor * m/2;
Vref_B_THI = (Vref_B_THI - offset) * scale_factor * m/2;
Vref_C_THI = (Vref_C_THI - offset) * scale_factor * m/2;

PWM_A_THI = double(Vref_A_THI > carrier);
PWM_B_THI = double(Vref_B_THI > carrier);
PWM_C_THI = double(Vref_C_THI > carrier);

%% Cálculo das Tensões de Saída e Correntes
fprintf('Calculando tensões e correntes de saída...\n');

% Tensões de fase para SPWM
Van_SPWM = (PWM_A_SPWM - 0.5) * Vdc;
Vbn_SPWM = (PWM_B_SPWM - 0.5) * Vdc;
Vcn_SPWM = (PWM_C_SPWM - 0.5) * Vdc;

% Tensões de linha
Vab_SPWM = Van_SPWM - Vbn_SPWM;
Vbc_SPWM = Vbn_SPWM - Vcn_SPWM;
Vca_SPWM = Vcn_SPWM - Van_SPWM;

% Simulação da carga RL trifásica
fprintf('Simulando carga RL trifásica...\n');
ia_SPWM = zeros(size(t));
ib_SPWM = zeros(size(t));
ic_SPWM = zeros(size(t));

for k = 2:length(t)
    dt = t(k) - t(k-1);
    % Fase A
    dia = (Van_SPWM(k) - R_load*ia_SPWM(k-1)) / L_load;
    ia_SPWM(k) = ia_SPWM(k-1) + dia * dt;
    
    % Fase B
    dib = (Vbn_SPWM(k) - R_load*ib_SPWM(k-1)) / L_load;
    ib_SPWM(k) = ib_SPWM(k-1) + dib * dt;
    
    % Fase C
    dic = (Vcn_SPWM(k) - R_load*ic_SPWM(k-1)) / L_load;
    ic_SPWM(k) = ic_SPWM(k-1) + dic * dt;
end


%% Análise de Harmônicas - FFT
fprintf('Realizando análise harmônica...\n');

% Janela para análise em regime permanente
t_steady = t(t > 0.04);  % Ignora transitório inicial
idx_steady = find(t > 0.04, 1);

% FFT da tensão de linha Vab
N_fft = 2^nextpow2(length(t_steady));
Y = fft(Vab_SPWM(idx_steady:end), N_fft);
P2 = abs(Y/N_fft);
P1 = P2(1:N_fft/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f_fft = f_chaveamento*(0:(N_fft/2))/N_fft;

% Encontrar harmônicas significativas
fundamental_idx = find(f_fft >= f_fundamental-1 & f_fft <= f_fundamental+1, 1);
fundamental_mag = P1(fundamental_idx);

% THD cálculo
harmonics_power = sum(P1(2:end).^2) - fundamental_mag^2; % Remove fundamental
THD = sqrt(harmonics_power) / fundamental_mag * 100;

fprintf('Resultados da Análise Harmônica:\n');
fprintf('  Tensão Fundamental: %.2f V\n', fundamental_mag);
fprintf('  THD: %.2f%%\n', THD);

%% Visualização dos Resultados
fprintf('Gerando gráficos...\n');

figure('Position', [100, 100, 1400, 1000]);

% 1. Sinais de Referência e Portadora
subplot(4,3,1);
plot(t(1:2000), Vref_A(1:2000), 'r', 'LineWidth', 1.5); hold on;
plot(t(1:2000), Vref_B(1:2000), 'g', 'LineWidth', 1.5);
plot(t(1:2000), Vref_C(1:2000), 'b', 'LineWidth', 1.5);
plot(t(1:2000), carrier(1:2000), 'k--', 'LineWidth', 1);
title('Sinais de Referência e Portadora');
xlabel('Tempo [s]'); ylabel('Amplitude [pu]');
legend('V_{ref,A}', 'V_{ref,B}', 'V_{ref,C}', 'Portadora', 'Location', 'best');
grid on;

% 2. Sinais PWM Fase A - Comparação de Técnicas
subplot(4,3,2);
plot(t(1:1000), PWM_A_SPWM(1:1000), 'r', 'LineWidth', 1); hold on;
plot(t(1:1000), PWM_A_SVPWM(1:1000), 'b', 'LineWidth', 1);
plot(t(1:1000), PWM_A_THI(1:1000), 'g', 'LineWidth', 1);
title('Sinais PWM Fase A - Técnicas Comparadas');
xlabel('Tempo [s]'); ylabel('Estado');
legend('SPWM', 'SVPWM', 'THIPWM', 'Location', 'best');
grid on;
ylim([-0.1 1.1]);

% 3. Tensão de Linha Vab
subplot(4,3,3);
plot(t(1:2000), Vab_SPWM(1:2000), 'b', 'LineWidth', 1.5);
title('Tensão de Linha V_{ab}');
xlabel('Tempo [s]'); ylabel('Tensão [V]');
grid on;

% 4. Correntes de Fase
subplot(4,3,4);
plot(t, ia_SPWM, 'r', 'LineWidth', 1.5); hold on;
plot(t, ib_SPWM, 'g', 'LineWidth', 1.5);
plot(t, ic_SPWM, 'b', 'LineWidth', 1.5);
title('Correntes de Fase');
xlabel('Tempo [s]'); ylabel('Corrente [A]');
legend('i_a', 'i_b', 'i_c', 'Location', 'best');
grid on;

% 5. Espectro Harmônico
subplot(4,3,5);
stem(f_fft(1:200), P1(1:200), 'filled', 'LineWidth', 1.5);
title('Espectro Harmônico - Tensão V_{ab}');
xlabel('Frequência [Hz]'); ylabel('Amplitude [V]');
xlim([0 2000]);
grid on;

% 6. Diagrama Vetorial
subplot(4,3,6);
theta = 0:pi/30:2*pi;
plot(cos(theta), sin(theta), 'k--'); hold on;
plot(V_alpha(1:100:end), V_beta(1:100:end), 'ro', 'MarkerSize', 3);
quiver(0, 0, mean(V_alpha(1000:2000)), mean(V_beta(1000:2000)), 'b', 'LineWidth', 2);
title('Diagrama Vetorial αβ');
xlabel('α'); ylabel('β');
axis equal; grid on;

% 7. THD vs Índice de Modulação
subplot(4,3,7);
m_range = 0.1:0.1:1.2;
THD_range = zeros(size(m_range));
for i = 1:length(m_range)
    Vref_test = m_range(i) * sin(2*pi*f_ref*t);
    PWM_test = double(Vref_test > carrier);
    Vout_test = (PWM_test - 0.5) * Vdc;
    Y_test = fft(Vout_test(idx_steady:end), N_fft);
    P2_test = abs(Y_test/N_fft);
    P1_test = P2_test(1:N_fft/2+1);
    P1_test(2:end-1) = 2*P1_test(2:end-1);
    fund_idx = find(f_fft >= f_fundamental-1 & f_fft <= f_fundamental+1, 1);
    fund_mag = P1_test(fund_idx);
    harmonics_power = sum(P1_test(2:end).^2) - fund_mag^2;
    THD_range(i) = sqrt(harmonics_power) / fund_mag * 100;
end
plot(m_range, THD_range, 'b-o', 'LineWidth', 2);
title('THD vs Índice de Modulação');
xlabel('Índice de Modulação m'); ylabel('THD [%]');
grid on;

% 8. Eficiência das Técnicas
subplot(4,3,8);
techniques = {'SPWM', 'SVPWM', 'THIPWM'};
fundamental_utilization = [0.5, 0.577, 0.577]; % Utilização da tensão DC
bar(fundamental_utilization, 'FaceColor', [0.2 0.6 0.8]);
set(gca, 'XTickLabel', techniques);
title('Utilização da Tensão DC');
ylabel('V_{fund} / V_{dc}');
grid on;

% 9. Formas de Onda Completas
subplot(4,3,9);
t_zoom = t(t > 0.05 & t < 0.06);
idx_zoom = find(t > 0.05 & t < 0.06);
plot(t_zoom, Vab_SPWM(idx_zoom), 'b', 'LineWidth', 1.5); hold on;
plot(t_zoom, ia_SPWM(idx_zoom)*50, 'r', 'LineWidth', 1.5); % Escalada para visualização
title('Tensão e Corrente (Zoom)');
xlabel('Tempo [s]'); ylabel('Tensão [V] / Corrente [A×50]');
legend('V_{ab}', 'i_a × 50', 'Location', 'best');
grid on;
