%% MASTER GPC + PWM CONTROL FOR 3-PHASE VSI - DUAL αβ TRACKING FROM DQ
% ═══════════════════════════════════════════════════════════════════════
% CORRIGIDO: GPC agora rastreia AMBOS α e β componentes
% Dois controladores GPC independentes: um para α, outro para β
% ═══════════════════════════════════════════════════════════════════════
clear; clc; close all;

%% Step 1: Clark Transformation Matrices
K_clark = 2/3 * [1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2];
K_clark_inv = [1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2];
I_2 = eye(2);

%% Step 2: System Parameters
Vdc = 60;               % DC link voltage [V]
R = 3.5;                % Resistance [Ω]
L = 0.1;                % Inductance [H]

% LC Filter Parameters
Rf = 0.1;               % Filter resistance [Ω]
Lf = 1e-3;              % Filter inductance [H]
Cf = 10e-6;             % Filter capacitance [F]

% Time parameters
f_sw = 10e3;            % Switching frequency [Hz]
f_fund = 60;            % Fundamental frequency [Hz]
Ts = 1/f_sw;            % Sampling time [s]
Tsim = 0.1;             % Simulation time [s]

fprintf('╔═══════════════════════════════════════════════════════╗\n');
fprintf('║  GPC + PWM CONTROL WITH DUAL αβ TRACKING (FROM DQ)   ║\n');
fprintf('╚═══════════════════════════════════════════════════════╝\n\n');

fprintf('=== SYSTEM PARAMETERS ===\n');
fprintf('Vdc = %.0f V, R = %.1f Ω, L = %.4f H\n', Vdc, R, L);
fprintf('f_sw = %.1f kHz, f_fund = %.0f Hz\n', f_sw/1e3, f_fund);
fprintf('Ts = %.6f s, Tsim = %.2f s\n\n', Ts, Tsim);

%% Step 3: Time Vector
t = 0:Ts:Tsim;
nit = length(t);
nin = 20;

fprintf('Simulation: nit=%d samples (%.2f ms)\n\n', nit, Tsim*1000);

%% Step 4: Discrete State-Space Model
% Continuous model for LC filter
A_cont = [-Rf/Lf,  -1/Lf;
           1/Cf,    0   ];
B_cont = [1/Lf; 0];
C_cont = [1 0];
D_cont = 0;

% Discretize
sysc = ss(A_cont, B_cont, C_cont, D_cont);
sysd = c2d(sysc, Ts, 'zoh');

% Extract transfer function
Tz = tf(sysd);
A_tf = cell2mat(Tz.den);   % A(z^-1)
B_tf = cell2mat(Tz.num);   % B(z^-1)

fprintf('=== DISCRETE MODEL ===\n');
fprintf('A(z⁻¹) = [%.6f, %.6f]\n', A_tf);
fprintf('B(z⁻¹) = [%.6f, %.6f]\n\n', B_tf);

%% Step 5: GPC Parameters
delta = 1;              % Output weight
lambda = 0.01;          % Control weight

d = 0;                  % Delay
N1 = 1;                 % Min horizon
N2 = 10;                % Max horizon
Nu = 3;                 % Control horizon
N = N2 - N1 + 1;        % Prediction length

na = length(A_tf) - 1;
nb = length(B_tf) - 1;

Bq = [0.001];           % Disturbance model
dq = 0;

fprintf('=== GPC CONFIGURATION ===\n');
fprintf('Horizons: N1=%d, N2=%d, Nu=%d → N=%d\n', N1, N2, Nu, N);
fprintf('Model orders: na=%d, nb=%d, d=%d\n', na, nb, d);
fprintf('Weights: δ=%.1f, λ=%.3f\n\n', delta, lambda);

%% Step 6: Build GPC Matrices
% Incorporate delay
Btil = conv(B_tf, [zeros(1,d) 1]);
Bqtil = conv(Bq, [zeros(1,dq) 1]);

% Diophantine equation
[E, F_poly] = diofantina(conv(A_tf, [1 -1]), N1, N2);

% Dynamic matrices
G = zeros(N, Nu);
H = zeros(N, nb + d);
Gq = zeros(N, 0);  % Not used
Hq = zeros(N, 0);  % Not used

for i = N1:N2
    EjB = conv(E(i-N1+1, 1:i), Btil);
    
    G(i-N1+1, 1:min(i, Nu)) = EjB(i:-1:max(1, i-Nu+1));
    
    if i < length(EjB)
        H(i-N1+1, 1:min(nb+d, length(EjB)-i)) = EjB(i+1:min(end, i+nb+d));
    end
end

G = G(:, 1:Nu);

fprintf('=== GPC MATRICES ===\n');
fprintf('G: %dx%d, H: %dx%d, F: %dx%d\n', size(G), size(H), size(F_poly));

%% Step 7: GPC Gain Calculation
Qu = lambda * eye(Nu);
Qe = delta * eye(N);

K_gpc = (G'*Qe*G + Qu) \ (G'*Qe);
K_gpc1 = K_gpc(1, :);  % First row (receding horizon)

fprintf('K_gpc: %dx%d\n', size(K_gpc));
fprintf('K_gpc1 (used): %dx%d\n\n', size(K_gpc1));

%% Step 8: Initialize Variables
% ═══════════════════════════════════════════════════════════════════════
% DUAL TRACKING: Separate variables for α and β components
% ═══════════════════════════════════════════════════════════════════════

% α-axis variables
saidas_alpha = zeros(nit, 1);       % α output
entradas_alpha = zeros(nit, 1);     % α control
du_alpha = zeros(nit, 1);           % α increments
refs_alpha = zeros(nit, 1);         % α reference

% β-axis variables
saidas_beta = zeros(nit, 1);        % β output
entradas_beta = zeros(nit, 1);      % β control
du_beta = zeros(nit, 1);            % β increments
refs_beta = zeros(nit, 1);          % β reference

% DQ references
Vd_ref = zeros(nit, 1);
Vq_ref = zeros(nit, 1);

% αβ references (from DQ)
V_alpha_ref = zeros(nit, 1);
V_beta_ref = zeros(nit, 1);

% PWM signals
pwm_a = zeros(nit, 1);
pwm_b = zeros(nit, 1);
pwm_c = zeros(nit, 1);

% Carrier
carrier = sawtooth(2*pi*f_sw*t, 0.5);

%% Step 9: Generate DQ References
fprintf('Generating DQ references...\n');

V_ref_peak = 5;         % Peak voltage [V]
f_ref = 60;             % Reference frequency [Hz]

for k = 1:nit
    if k > nin
        % Vd = constant (DC), Vq = 0 for simple test
        Vd_ref(k) = V_ref_peak;
        Vq_ref(k) = 0;
        
        % For rotating reference:
        % Vd_ref(k) = V_ref_peak * cos(2*pi*10*t(k));
        % Vq_ref(k) = V_ref_peak * sin(2*pi*10*t(k));
    end
end

% Inverse Park transformation: DQ → αβ
for k = 1:nit
    theta = 2*pi*f_ref*t(k);
    
    % [α]   [cos(θ)  -sin(θ)] [d]
    % [β] = [sin(θ)   cos(θ)] [q]
    V_alpha_ref(k) = Vd_ref(k)*cos(theta) - Vq_ref(k)*sin(theta);
    V_beta_ref(k)  = Vd_ref(k)*sin(theta) + Vq_ref(k)*cos(theta);
    
    % Set GPC references
    refs_alpha(k) = V_alpha_ref(k);
    refs_beta(k)  = V_beta_ref(k);
end

fprintf('DQ → αβ transformation complete.\n');
fprintf('  Max V_α: %.2f V, Max V_β: %.2f V\n\n', ...
        max(abs(V_alpha_ref)), max(abs(V_beta_ref)));

%% Step 10: Main Simulation Loop with DUAL GPC
fprintf('Starting DUAL GPC + PWM simulation...\n');

for k = nin:nit-1
    
    %% ═══════════════════════════════════════════════════════════════════
    %  PLANTA: Simula ambos os canais α e β independentemente
    %% ═══════════════════════════════════════════════════════════════════
    
    if k > na && k > nb+d
        % α-axis plant
        past_outputs_alpha = saidas_alpha(k:-1:k-na);
        past_inputs_alpha = entradas_alpha(k-d:-1:k-nb-d);
        saidas_alpha(k+1) = -A_tf(2:end) * past_outputs_alpha(2:end) + ...
                             B_tf * past_inputs_alpha;
        
        % β-axis plant
        past_outputs_beta = saidas_beta(k:-1:k-na);
        past_inputs_beta = entradas_beta(k-d:-1:k-nb-d);
        saidas_beta(k+1) = -A_tf(2:end) * past_outputs_beta(2:end) + ...
                            B_tf * past_inputs_beta;
    else
        saidas_alpha(k+1) = saidas_alpha(k);
        saidas_beta(k+1) = saidas_beta(k);
    end
    
    % Measurement noise
    noise_alpha = 0.01 * randn;
    noise_beta = 0.01 * randn;
    
    %% ═══════════════════════════════════════════════════════════════════
    %  GPC CONTROLADOR ALPHA (Eixo α)
    %% ═══════════════════════════════════════════════════════════════════
    
    if k > N2
        % Reference vector
        R_alpha = refs_alpha(k) * ones(N, 1);
        
        % Free response
        f_alpha = zeros(N, 1);
        for j = 1:N
            if j <= size(F_poly, 1)
                f_alpha(j) = F_poly(j, :) * saidas_alpha(k:-1:max(1, k-na));
            else
                f_alpha(j) = saidas_alpha(k);
            end
        end
        
        % Add past control increments
        if ~isempty(H) && k > nb+d
            past_du_alpha = du_alpha(k-1:-1:max(1, k-nb-d));
            if length(past_du_alpha) < size(H, 2)
                past_du_alpha = [past_du_alpha; zeros(size(H,2)-length(past_du_alpha), 1)];
            end
            for j = 1:N
                f_alpha(j) = f_alpha(j) + H(j, :) * past_du_alpha;
            end
        end
        
        % Control law
        du_alpha(k) = K_gpc1 * (R_alpha - f_alpha);
        entradas_alpha(k) = entradas_alpha(k-1) + du_alpha(k);
        
        % Saturation
        max_v = Vdc / sqrt(3);
        entradas_alpha(k) = max(min(entradas_alpha(k), max_v), -max_v);
    else
        entradas_alpha(k) = refs_alpha(k);
        du_alpha(k) = 0;
    end
    
    %% ═══════════════════════════════════════════════════════════════════
    %  GPC CONTROLADOR BETA (Eixo β)
    %% ═══════════════════════════════════════════════════════════════════
    
    if k > N2
        % Reference vector
        R_beta = refs_beta(k) * ones(N, 1);
        
        % Free response
        f_beta = zeros(N, 1);
        for j = 1:N
            if j <= size(F_poly, 1)
                f_beta(j) = F_poly(j, :) * saidas_beta(k:-1:max(1, k-na));
            else
                f_beta(j) = saidas_beta(k);
            end
        end
        
        % Add past control increments
        if ~isempty(H) && k > nb+d
            past_du_beta = du_beta(k-1:-1:max(1, k-nb-d));
            if length(past_du_beta) < size(H, 2)
                past_du_beta = [past_du_beta; zeros(size(H,2)-length(past_du_beta), 1)];
            end
            for j = 1:N
                f_beta(j) = f_beta(j) + H(j, :) * past_du_beta;
            end
        end
        
        % Control law
        du_beta(k) = K_gpc1 * (R_beta - f_beta);
        entradas_beta(k) = entradas_beta(k-1) + du_beta(k);
        
        % Saturation
        max_v = Vdc / sqrt(3);
        entradas_beta(k) = max(min(entradas_beta(k), max_v), -max_v);
    else
        entradas_beta(k) = refs_beta(k);
        du_beta(k) = 0;
    end
    
    %% ═══════════════════════════════════════════════════════════════════
    %  SPACE VECTOR MODULATION (SVM)
    %  Usa AMBOS V_α e V_β dos controladores GPC
    %% ═══════════════════════════════════════════════════════════════════
    
    % Voltage vector in complex form
    V_ref_complex = entradas_alpha(k) + 1j*entradas_beta(k);
    
    % Magnitude and angle
    V_mag = abs(V_ref_complex);
    theta_svm = angle(V_ref_complex);
    theta_svm = mod(theta_svm, 2*pi);
    
    % Identify sector (1-6)
    sector = floor(theta_svm / (pi/3)) + 1;
    if sector > 6, sector = 6; end
    
    % Switching times
    alpha_angle = mod(theta_svm, pi/3);
    
    T1 = Ts * V_mag * sin(pi/3 - alpha_angle) / Vdc;
    T2 = Ts * V_mag * sin(alpha_angle) / Vdc;
    T0 = Ts - T1 - T2;
    
    % Handle overmodulation
    if T0 < 0
        T_total = T1 + T2;
        T1 = T1 * Ts / T_total;
        T2 = T2 * Ts / T_total;
        T0 = 0;
    end
    
    % Duty cycles based on sector
    switch sector
        case 1
            Ta = (T1 + T2 + T0/2) / Ts;
            Tb = (T2 + T0/2) / Ts;
            Tc = T0/2 / Ts;
        case 2
            Ta = (T1 + T0/2) / Ts;
            Tb = (T1 + T2 + T0/2) / Ts;
            Tc = T0/2 / Ts;
        case 3
            Ta = T0/2 / Ts;
            Tb = (T1 + T2 + T0/2) / Ts;
            Tc = (T2 + T0/2) / Ts;
        case 4
            Ta = T0/2 / Ts;
            Tb = (T1 + T0/2) / Ts;
            Tc = (T1 + T2 + T0/2) / Ts;
        case 5
            Ta = (T2 + T0/2) / Ts;
            Tb = T0/2 / Ts;
            Tc = (T1 + T2 + T0/2) / Ts;
        case 6
            Ta = (T1 + T2 + T0/2) / Ts;
            Tb = T0/2 / Ts;
            Tc = (T1 + T0/2) / Ts;
    end
    
    % Limit duty cycles
    Ta = max(min(Ta, 1), 0);
    Tb = max(min(Tb, 1), 0);
    Tc = max(min(Tc, 1), 0);
    
    % PWM signals
    pwm_a(k) = Ta;
    pwm_b(k) = Tb;
    pwm_c(k) = Tc;
end

fprintf('Simulation completed!\n\n');

%% Step 11: Performance Analysis
fprintf('═══════════════════════════════════════════════════════\n');
fprintf('           PERFORMANCE ANALYSIS                        \n');
fprintf('═══════════════════════════════════════════════════════\n\n');

steady_start = find(t > 0.03, 1);

% α-axis metrics
error_alpha = saidas_alpha - refs_alpha;
rms_error_alpha = rms(error_alpha(steady_start:end));
tracking_alpha = 100 * (1 - rms_error_alpha / rms(refs_alpha(steady_start:end)));

% β-axis metrics
error_beta = saidas_beta - refs_beta;
rms_error_beta = rms(error_beta(steady_start:end));
tracking_beta = 100 * (1 - rms_error_beta / rms(refs_beta(steady_start:end)));

% Combined metrics
total_error = sqrt(error_alpha.^2 + error_beta.^2);
rms_total = rms(total_error(steady_start:end));

fprintf('╔═══════════════════════════════════════════════════════╗\n');
fprintf('║              α-AXIS PERFORMANCE                       ║\n');
fprintf('╠═══════════════════════════════════════════════════════╣\n');
fprintf('║  RMS Error:          %.4f V                       ║\n', rms_error_alpha);
fprintf('║  Tracking Accuracy:  %.2f%%                       ║\n', tracking_alpha);
fprintf('║  Control Effort:     %.4f                         ║\n', rms(du_alpha(nin:end)));
fprintf('╚═══════════════════════════════════════════════════════╝\n\n');

fprintf('╔═══════════════════════════════════════════════════════╗\n');
fprintf('║              β-AXIS PERFORMANCE                       ║\n');
fprintf('╠═══════════════════════════════════════════════════════╣\n');
fprintf('║  RMS Error:          %.4f V                       ║\n', rms_error_beta);
fprintf('║  Tracking Accuracy:  %.2f%%                       ║\n', tracking_beta);
fprintf('║  Control Effort:     %.4f                         ║\n', rms(du_beta(nin:end)));
fprintf('╚═══════════════════════════════════════════════════════╝\n\n');

fprintf('╔═══════════════════════════════════════════════════════╗\n');
fprintf('║           COMBINED αβ PERFORMANCE                     ║\n');
fprintf('╠═══════════════════════════════════════════════════════╣\n');
fprintf('║  Total RMS Error:    %.4f V                       ║\n', rms_total);
fprintf('╚═══════════════════════════════════════════════════════╝\n\n');

%% Step 12: Comprehensive Plotting
figure('Position', [50, 50, 1600, 1000]);

time_idx = nin:nit;

% Subplot 1: α-axis tracking
subplot(3, 4, 1);
plot(t(time_idx), saidas_alpha(time_idx), 'b-', 'LineWidth', 2); hold on;
plot(t(time_idx), refs_alpha(time_idx), 'r--', 'LineWidth', 2);
title('α-Axis Voltage Tracking');
xlabel('Time [s]'); ylabel('V_α [V]');
legend('Output', 'Reference', 'Location', 'best');
grid on;

% Subplot 2: β-axis tracking
subplot(3, 4, 2);
plot(t(time_idx), saidas_beta(time_idx), 'b-', 'LineWidth', 2); hold on;
plot(t(time_idx), refs_beta(time_idx), 'r--', 'LineWidth', 2);
title('β-Axis Voltage Tracking');
xlabel('Time [s]'); ylabel('V_β [V]');
legend('Output', 'Reference', 'Location', 'best');
grid on;

% Subplot 3: DQ references
subplot(3, 4, 3);
plot(t, Vd_ref, 'b-', 'LineWidth', 2); hold on;
plot(t, Vq_ref, 'r--', 'LineWidth', 2);
title('DQ Reference Voltages');
xlabel('Time [s]'); ylabel('Voltage [V]');
legend('V_d', 'V_q', 'Location', 'best');
grid on;

% Subplot 4: αβ vector trajectory
subplot(3, 4, 4);
plot_idx = nin:5:nit;
plot(entradas_alpha(plot_idx), entradas_beta(plot_idx), 'bo', 'MarkerSize', 3); hold on;
plot(refs_alpha(plot_idx), refs_beta(plot_idx), 'rx', 'MarkerSize', 4);
theta_c = 0:0.01:2*pi;
plot(V_ref_peak*cos(theta_c), V_ref_peak*sin(theta_c), 'k--');
title('αβ Vector Trajectory');
xlabel('V_α [V]'); ylabel('V_β [V]');
legend('GPC Output', 'Reference', 'Max Circle', 'Location', 'best');
axis equal; grid on;

% Subplot 5: α control input
subplot(3, 4, 5);
plot(t(time_idx), entradas_alpha(time_idx), 'g-', 'LineWidth', 2);
title('α-Axis Control Input');
xlabel('Time [s]'); ylabel('u_α [V]');
grid on;

% Subplot 6: β control input
subplot(3, 4, 6);
plot(t(time_idx), entradas_beta(time_idx), 'm-', 'LineWidth', 2);
title('β-Axis Control Input');
xlabel('Time [s]'); ylabel('u_β [V]');
grid on;

% Subplot 7: α error
subplot(3, 4, 7);
plot(t(time_idx), error_alpha(time_idx), 'r-', 'LineWidth', 2);
title('α-Axis Tracking Error');
xlabel('Time [s]'); ylabel('Error [V]');
grid on;

% Subplot 8: β error
subplot(3, 4, 8);
plot(t(time_idx), error_beta(time_idx), 'r-', 'LineWidth', 2);
title('β-Axis Tracking Error');
xlabel('Time [s]'); ylabel('Error [V]');
grid on;

% Subplot 9: PWM signals
subplot(3, 4, 9);
zoom_start = round(0.05/Ts);
zoom_end = round(0.052/Ts);
t_zoom = t(zoom_start:zoom_end);
plot(t_zoom, pwm_a(zoom_start:zoom_end), 'r', 'LineWidth', 2); hold on;
plot(t_zoom, pwm_b(zoom_start:zoom_end), 'g', 'LineWidth', 2);
plot(t_zoom, pwm_c(zoom_start:zoom_end), 'b', 'LineWidth', 2);
title('PWM Signals (50-52 ms)');
xlabel('Time [s]'); ylabel('Duty Cycle');
legend('Phase A', 'Phase B', 'Phase C', 'Location', 'best');
grid on; ylim([-0.1, 1.1]);

% Subplot 10: Control increments
subplot(3, 4, 10);
plot(t(time_idx), du_alpha(time_idx), 'b-', 'LineWidth', 1.5); hold on;
plot(t(time_idx), du_beta(time_idx), 'r-', 'LineWidth', 1.5);
title('Control Increments');
xlabel('Time [s]'); ylabel('\Delta u');
legend('\Delta u_α', '\Delta u_β', 'Location', 'best');
grid on;

% Subplot 11: Total error
subplot(3, 4, 11);
plot(t(time_idx), total_error(time_idx), 'k-', 'LineWidth', 2);
title('Total Vector Error |\vec{e}|');
xlabel('Time [s]'); ylabel('Error [V]');
grid on;

% Subplot 12: Performance comparison
subplot(3, 4, 12);
metrics = [rms_error_alpha, rms_error_beta, tracking_alpha, tracking_beta];
bar_labels = {'α RMS', 'β RMS', 'α Track%', 'β Track%'};
bar(metrics);
set(gca, 'XTickLabel', bar_labels);
title('Performance Metrics');
ylabel('Value');
grid on;

sgtitle('DUAL GPC + PWM: Independent αβ Tracking from DQ References', ...
        'FontSize', 16, 'FontWeight', 'bold');

fprintf('Plotting complete!\n');

%% Diophantine Function
function [E, F] = diofantina(A, N1, N2)
    nA = length(A) - 1;
    
    % Initialize F
    f = zeros(N2+1, nA);
    f(1, 1) = 1;
    
    % Calculate F polynomials
    for j = 1:N2
        for i = 1:nA-1
            f(j+1, i) = f(j, i+1) - f(j, 1) * A(i+1);
        end
        if nA > 0
            f(j+1, nA) = -f(j, 1) * A(nA+1);
        end
    end
    
    F = f(N1+1:N2+1, :);
    
    % Calculate E polynomials
    E = zeros(N2-N1+1, N2);
    e = zeros(1, N2);
    e(1) = 1;
    E(1, 1) = e(1);
    
    for i = 2:N2
        e(i) = f(i, 1);
        if i >= N1
            E(i-N1+1, 1:i) = e(1:i);
        end
    end
end
