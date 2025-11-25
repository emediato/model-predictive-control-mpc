%% MASTER GPC + PWM CONTROL FOR 3-PHASE VSI - DQ VOLTAGE REFERENCE
% Complete integration of GPC predictive control with PWM modulation
clear; clc; close all;

%% Step 1: Clark Transformation Matrices
K_clark = 2/3 * [1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2]; K=K_clark;
K_clark_inv = [1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2];
I_2 = eye(2);

%% Step 2: System Parameters Definition
% Hardware parameters
Vdc = 60;              % DC link voltage [V]
R = 3.5;                % Resistance [Ω]
L = 0.1;             % Inductance [H]

% LC Filter Parameters
Rf = 0.1;               % Filter resistance [Ω]
Lf = 1e-3;              % Filter inductance [H]
Cf = 10e-6;             % Filter capacitance [F]

% Time parameters
f_sw = 10e3;            % Switching frequency [Hz]
f_fund = 60;            % Fundamental frequency [Hz]
Ts = 1/f_sw;            % Sampling time [s]
Tsim = 0.1;             % Simulation time [s]

fprintf('=== SYSTEM PARAMETERS ===\n');
fprintf('Vdc = %.0f V, R = %.1f Ω, L = %.4f H\n', Vdc, R, L);
fprintf('f_sw = %.1f kHz, f_fund = %.0f Hz, Ts = %.6f s\n\n', f_sw/1e3, f_fund, Ts);

%% Step 3: Simulation Time Setup
t = 0:Ts:Tsim;
nit = length(t);        % Total number of iterations
nin = 20;               % Initial settling samples

fprintf('Simulation parameters: nit=%d, nin=%d\n\n', nit, nin);

% state matrix
F = - R/L * I_2;
G = (Vdc/(2*L)) * K;
A = expm(F*Ts); % matrix exponential
B = -F^(-1) * (I_2 - A) * G;

 
%% Step 4: Discrete State-Space Model for GPC
% Continuous-time matrices for LC filter in αβ frame
A_cont = [-Rf/Lf,  -1/Lf;
           1/Cf,    0   ];
B_cont = [1/Lf; 0];
C_cont = [1 0];
D_cont = 0;

% Discretize using exact method
sysc = ss(A_cont, B_cont, C_cont, D_cont);
sysd = c2d(sysc, Ts, 'zoh');
Ad = sysd.A;
Bd = sysd.B;

% Extract discrete transfer function for GPC
Tz = tf(sysd);
A_tf = cell2mat(Tz.den);   % A(z^-1) coefficients
B_tf = cell2mat(Tz.num);   % B(z^-1) coefficients
B = B_tf;
A = A_tf;

Bq = [0.001]; % numerador do modelo da perturbação do processo


fprintf('Discrete Model Coefficients:\n');
fprintf('A = [%.6f, %.6f]\n', A_tf);
fprintf('B = [%.6f, %.6f]\n\n', B_tf);


delta = 1;              % Output error weight

% ponderação do esforço de controle
lambda = 0.01; % switch
lambda_i = 10 ; % current

na = length(A_tf) - 1;
nb = length(B_tf) - 1;

fprintf('=== GPC PARAMETERS ===\n');

%% Step 5: GPC Parameters Setup
d = 0;                  % Input delay
N1 = 1;                 % Minimum prediction horizon
N2=10;
N2 = N2-N1+1;                % Maximum prediction horizon  
Nu = 3;                 % Control horizon
N = N2 - N1 + 1;        % Prediction horizon length
Nq = 0;%  horizonte da ação antecipativa


fprintf('N1=%d, N2=%d, Nu=%d, N=%d\n', N1, N2, Nu, N);
fprintf('na=%d, nb=%d, d=%d\n', na, nb, d);
fprintf('delta=%.1f, lambda=%.1f\n\n', delta, lambda);

%% Build GPC Matrices
% Incorporate delay into B polynomial

d=0; % atraso em relação à manipulada
dq = 0; % atraso em relação à perturbação


d=0; % atraso em relação à manipulada
dq = 0; % atraso em relação à perturbação

% montagem das matrizes
Btil = conv(B,[zeros(1,d) 1]); % incorporação do atraso no numerador B
Bqtil = conv(Bq,[zeros(1,dq) 1]); % incorporação do atraso no numerador Bq

% Calculate Diophantine polynomials
[E, F_poly] = diofantina(conv(A_tf, [1 -1]), N1, N2);

% Build dynamic matrix G and free response matrix H
G = zeros(N, Nu);
H = zeros(N, nb + d);


[E,F] = diofantina(conv2(A,[1 -1]),N1,N2); % obtenção dos polinômios Ej, Fj

for i=N1:N2
    EjB = conv2(E(i-N1+1,1:i),Btil);
    EjBq = conv2(E(i-N1+1,1:i),Bqtil);
    
    G(i-N1+1,1:i) = EjB(i:-1:1);    
    H(i-N1+1,:) = EjB(i+1:end);

    Gq(i-N1+1,1:i) = EjBq(i:-1:1);    
    Hq(i-N1+1,:) = EjBq(i+1:end);
    
end
G = G(:,1:Nu);
Gq = Gq(:,1:Nq);


G,F,H,Gq,Hq

%% obtenção do ganho irrestrito
Qu = lambda*eye(Nu); % matriz de ponderação dos incrementos de controle
Qe = delta*eye(N);  % matriz de ponderaão dos erros futuros

Kmpc = (G'*Qe*G+Qu)\G'*Qe;
Kmpc1 = Kmpc(1,:);

%% inicialização vetores
nin = 5;
nit = 100 + nin; % número de iterações da simulação

K_gpc = Kmpc;
K_gpc1 = K_gpc(1, :);  % First row for current control

fprintf('GPC Matrix Dimensions:\n');
fprintf('G: %dx%d, H: %dx%d\n', size(G,1), size(G,2), size(H,1), size(H,2));
fprintf('K_gpc: %dx%d\n\n', size(K_gpc,1), size(K_gpc,2));

%%  Initialize Simulation Variables
% System variables
saidas = zeros(nit, 1);     % System outputs (currents)
entradas = zeros(nit, 1);   % Control inputs (voltages)  
du = zeros(nit, 1);         % Control increments
refs = zeros(nit, 1);       % References

% DQ reference voltages
Vd_ref = zeros(nit, 1);     % D-axis voltage reference
Vq_ref = zeros(nit, 1);     % Q-axis voltage reference

% αβ reference voltages (from DQ transformation)
V_alpha_ref = zeros(nit, 1);
V_beta_ref = zeros(nit, 1);

% PWM signals
pwm_a = zeros(nit, 1);
pwm_b = zeros(nit, 1);
pwm_c = zeros(nit, 1);

% Current measurements
I_alpha_meas = zeros(nit, 1);
I_beta_meas = zeros(nit, 1);

% Carrier signal for PWM
carrier = sawtooth(2*pi*f_sw*t, 0.5);

%% Step 8: Reference Generation - DQ Voltage References
fprintf('Generating DQ voltage references...\n');

% Reference parameters
V_ref_peak = 5;       % Peak voltage reference [V]
f_ref = 60;             % Reference frequency [Hz]

for k = 1:nit
    if k > nin
        % Step reference in D-axis, zero Q-axis for simplicity
        Vd_ref(k) = V_ref_peak;
        Vq_ref(k) = 0;
        
        % For sinusoidal reference tracking:
        % Vd_ref(k) = V_ref_peak * cos(2*pi*f_ref*t(k));
        % Vq_ref(k) = V_ref_peak * sin(2*pi*f_ref*t(k));
    end
end

% Transform DQ to αβ frame (inverse Park transformation)
for k = 1:nit
    theta = 2*pi*f_ref*t(k);  % Angle for transformation
    
    % Inverse Park transformation: αβ to DQ
    V_alpha_ref(k) = Vd_ref(k)*cos(theta) - Vq_ref(k)*sin(theta);
    V_beta_ref(k) = Vd_ref(k)*sin(theta) + Vq_ref(k)*cos(theta);
    
    % Set current reference for GPC (track voltage)
    refs(k) = V_alpha_ref(k);  % GPC tracks α-component
end

%% Step 9: Perturbation Signal
perts = zeros(nit, 1);
perts(nin+50:end) = 0.5;  % Step perturbation

%% Step 10: Main GPC + PWM Simulation Loop
fprintf('Starting GPC + PWM simulation with DQ references...\n');

for k = nin:nit-1
    %% Step 10.1: System Model (Plant Simulation)
    % Simulate the discrete system using transfer function model
    if k > na && k > nb+d
        past_outputs = saidas(k:-1:k-na);
        past_inputs = entradas(k-d:-1:k-nb-d);
        
        saidas(k+1) = -A_tf(2:end) * past_outputs(2:end) + ...
                       B_tf * past_inputs;
    else
        saidas(k+1) = saidas(k);  % Maintain previous value
    end
    
    % Add measurement noise to simulate real sensors
    measurement_noise = 0.01 * randn;
    I_alpha_meas(k) = saidas(k) + measurement_noise;
    
    %% GPC Control Law
    if k > N2
        % Reference vector for prediction horizon
        R = refs(k) * ones(N, 1);
        
        % Free response calculation
        f = zeros(N, 1);
        current_output = I_alpha_meas(k);
        
        % Simplified free response prediction
        for j = 1:N
            f(j) = current_output;  % Basic constant prediction
        end
        
        % Add contribution from past control actions
        if ~isempty(H) && k > nb+d
            past_du = du(k-1:-1:max(1,k-nb-d));
            if length(past_du) < size(H,2)
                past_du = [past_du; zeros(size(H,2)-length(past_du), 1)];
            end
            f = f + H * past_du;
        end
        
        % Control increment calculation
        du(k) = K_gpc1 * (R - f);
        entradas(k) = entradas(k-1) + du(k);
        
        % Anti-windup saturation
        max_voltage = Vdc/sqrt(3);
        entradas(k) = max(min(entradas(k), max_voltage), -max_voltage);
    else
        % Initialization phase
        entradas(k) = refs(k);
        du(k) = 0;
    end
    
    %% Step 10.3: Space Vector Modulation (SVM)
    % Use GPC output as voltage reference for SVM
    V_ref_complex = entradas(k) + 1j * 0;  % β-component is zero for simplicity
    
    % Voltage magnitude and angle
    V_mag = abs(V_ref_complex);
    theta_svm = angle(V_ref_complex);
    theta_svm = mod(theta_svm, 2*pi);  % Normalize angle
    
    % Identify sector (1 to 6)
    sector = floor(theta_svm / (pi/3)) + 1;
    
    % Calculate switching times
    alpha = mod(theta_svm, pi/3);
    
    T1 = Ts * V_mag * sin(pi/3 - alpha) / Vdc;
    T2 = Ts * V_mag * sin(alpha) / Vdc;
    T0 = Ts - T1 - T2;
    
    % Calculate duty cycles based on sector
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
    
    % Ensure duty cycles are within [0,1]
    Ta = max(min(Ta, 1), 0);
    Tb = max(min(Tb, 1), 0);
    Tc = max(min(Tc, 1), 0);
    
    % Generate PWM signals by comparing with carrier
    pwm_a(k) = double(Ta > (carrier(k) + 1)/2);
    pwm_b(k) = double(Tb > (carrier(k) + 1)/2);
    pwm_c(k) = double(Tc > (carrier(k) + 1)/2);
end

fprintf('Simulation completed!\n');

%% Step 11: Performance Analysis
fprintf('\n=== PERFORMANCE ANALYSIS ===\n');

% Calculate performance metrics
steady_start = find(t > 0.02, 1);
steady_state_error = rms(saidas(steady_start:end) - refs(steady_start:end));
control_effort = rms(du(nin:end));
tracking_accuracy = 100 * (1 - steady_state_error / rms(refs(steady_start:end)));

fprintf('Steady-state RMS error: %.4f V\n', steady_state_error);
fprintf('Control effort (RMS): %.4f\n', control_effort);
fprintf('Tracking accuracy: %.2f%%\n', tracking_accuracy);

%% Step 12: Comprehensive Results Plotting
figure('Position', [50, 50, 1400, 1000]);
time_plot = 1:nit;
% Subplot 1: Voltage tracking performance
subplot(3, 3, 1);
time_plot_adjs = 80:nit;

plot(time_plot_adjs, saidas(time_plot_adjs), 'b-', 'LineWidth', 2); hold on;
plot(time_plot_adjs, refs(time_plot_adjs), 'r--', 'LineWidth', 2);
plot(time_plot_adjs, V_alpha_ref(time_plot_adjs), 'g-.', 'LineWidth', 1.5);
title('GPC Voltage Tracking Performance');
xlabel('Time [s]'); ylabel('Voltage [V]');
legend('System Output', 'GPC Reference', 'DQ α-Reference', 'Location', 'best');
grid on;

%%
% Subplot 2: DQ reference voltages
subplot(3, 3, 2);
plot(time_plot, Vd_ref, 'b-', 'LineWidth', 2); hold on;
plot(time_plot, Vq_ref, 'r--', 'LineWidth', 2);
title('DQ Reference Voltages');
xlabel('Time [s]'); ylabel('Voltage [V]');
legend('V_d Reference', 'V_q Reference', 'Location', 'best');
grid on;

% Subplot 3: Control input voltage
subplot(3, 3, 3);
plot(time_plot, entradas, 'g-', 'LineWidth', 2);
title('GPC Control Input Voltage');
xlabel('Time [s]'); ylabel('Voltage [V]');
grid on;

% Subplot 4: Control increments
subplot(3, 3, 4);
plot(time_plot, du, 'm-', 'LineWidth', 2);
title('Control Increments \Delta u');
xlabel('Time [s]'); ylabel('\Delta u');
grid on;

% Subplot 5: PWM signals (zoomed view)
subplot(3, 3, 5);
zoom_start = find(t > 0.05, 1);
zoom_end = find(t > 0.051, 1);

t_zoom = t(zoom_start:zoom_end);

plot(time_plot, pwm_a, 'r', 'LineWidth', 2); hold on;
plot(time_plot, pwm_b , 'g', 'LineWidth', 2);
plot(time_plot, pwm_c , 'b', 'LineWidth', 2);
% plot(time_plot, carrier , 'k--', 'LineWidth', 1);
title('PWM Signals (50-51 ms)');
xlabel('Time [s]'); ylabel('State');
legend('PWM A', 'PWM B', 'PWM C', 'Carrier', 'Location', 'best');
grid on;
ylim([-0.1 1.1]);

% Subplot 6: Tracking error
subplot(3, 3, 6);
error = saidas - refs;
plot(time_plot, error, 'r-', 'LineWidth', 2);
title('Tracking Error');
xlabel('Time [s]'); ylabel('Error [V]');
grid on;

% Subplot 7: αβ reference voltages
subplot(3, 3, 7);
plot(time_plot, V_alpha_ref, 'b-', 'LineWidth', 2); hold on;
plot(time_plot, V_beta_ref, 'r--', 'LineWidth', 2);
title('αβ Reference Voltages from DQ');
xlabel('Time [s]'); ylabel('Voltage [V]');
legend('V_{α,ref}', 'V_{β,ref}', 'Location', 'best');
grid on;

% Subplot 8: Performance metrics
subplot(3, 3, 8);
metrics = [steady_state_error, control_effort, max(abs(du)), tracking_accuracy];
metric_names = {'RMS Error', 'Control Effort', 'Max \Delta u', 'Accuracy %'};
bar(metrics);
set(gca, 'XTickLabel', metric_names);
title('Performance Metrics');
ylabel('Value');
grid on;

% Subplot 9: Voltage vector trajectory
subplot(3, 3, 9);
plot_idx = nin:10:nit;  % Plot every 10th sample for clarity
plot(entradas(plot_idx), zeros(size(plot_idx)), 'bo', 'MarkerSize', 3);
hold on;
plot(V_alpha_ref(plot_idx), V_beta_ref(plot_idx), 'rx', 'MarkerSize', 4);
theta_circle = 0:0.01:2*pi;
plot(V_ref_peak*cos(theta_circle), V_ref_peak*sin(theta_circle), 'k--');
title('Voltage Vector Trajectory');
xlabel('V_α [V]'); ylabel('V_β [V]');
legend('GPC Output', 'DQ Reference', 'Max Circle', 'Location', 'best');
axis equal;
grid on;

sgtitle('MASTER GPC + PWM CONTROL WITH DQ REFERENCES - 3-PHASE VSI', ...
        'FontSize', 16, 'FontWeight', 'bold');

%% Diofantina Function (Required for GPC)
function [E, F] = diofantina(A, N1, N2)
    % Diophantine equation solver for GPC
    % Solves: 1 = E_j(z^{-1})A(z^{-1}) + z^{-j}F_j(z^{-1})
    
    nA = length(A) - 1;
    
    % Initialize F matrix
    f = zeros(N2+1, nA);
    f(1, 1) = 1;  % F_0 = 1
    
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