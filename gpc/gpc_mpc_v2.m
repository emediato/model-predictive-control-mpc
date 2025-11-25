%% LC Filter State-Space Model for GPC with PWM Integration
% This script implements GPC control with PWM modulation for a 3-phase VSI
clear; clc; close all;

%% Step 1: Clark Transformation Matrix
K = 2/3 * [1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2];
K_inv = [1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2]; % Inverse Clark

I_2 = eye(2);

%% Step 2: System Parameters
% Hardware parameters
Vdc = 400;              % DC link voltage [V]
R = 3.5;                % Resistance [Ω]
L = 0.0025;             % Inductance [H]

% LC Filter Parameters (for more accurate model)
Rf = 0.1;               % Filter resistance [Ω]
Lf = 1e-3;              % Filter inductance [H]
Cf = 10e-6;             % Filter capacitance [F]

% Sampling and switching frequencies
f_sw = 10e3;            % Switching frequency [Hz]
f_fund = 60;            % Fundamental frequency [Hz]
Ts = 1/f_sw;            % Sampling time [s]
Tsim = 0.1;             % Simulation time [s]

fprintf('=== SYSTEM PARAMETERS ===\n');
fprintf('Vdc = %.0f V, R = %.1f Ω, L = %.4f H\n', Vdc, R, L);
fprintf('f_sw = %.1f kHz, f_fund = %.0f Hz, Ts = %.6f s\n\n', f_sw/1e3, f_fund, Ts);

%% Step 3: Discrete State-Space Model for GPC
% Continuous-time matrices for LC filter
A1 = [-Rf/Lf,  -1/Lf;
       1/Cf,    0   ];
B1 = [1/Lf; 0];
C1 = [1 0];
D1 = 0;

% Discretize using exact method
sysc = ss(A1, B1, C1, D1);
sysd = c2d(sysc, Ts, 'zoh');
Ad = sysd.A;
Bd = sysd.B;

% Extract discrete transfer function for GPC
Tz = tf(sysd);
A = cell2mat(Tz.den);   % A(z^-1) coefficients
B = cell2mat(Tz.num);   % B(z^-1) coefficients

fprintf('Discrete Model:\n');
fprintf('A = [%.6f, %.6f]\n', A);
fprintf('B = [%.6f, %.6f]\n\n', B);

%% Step 4: GPC Parameters
d = 0;                  % Input delay
N1 = 1;                 % Minimum prediction horizon
N2 = 10;                % Maximum prediction horizon  
Nu = 3;                 % Control horizon
N = N2 - N1 + 1;        % Prediction horizon length

delta = 1;              % Output error weight
lambda = 0.1;           % Control effort weight

fprintf('=== GPC PARAMETERS ===\n');
fprintf('N1=%d, N2=%d, Nu=%d, N=%d\n', N1, N2, Nu, N);
fprintf('delta=%.1f, lambda=%.1f\n\n', delta, lambda);

%% Step 5: Build GPC Matrices
na = length(A) - 1;
nb = length(B) - 1;

% Incorporate delay
Btil = conv(B, [zeros(1,d) 1]);

% Calculate Diophantine polynomials
[E, F] = diofantina(conv(A, [1 -1]), N1, N2);

% Build dynamic matrix G
G = zeros(N, Nu);
H = zeros(N, nb + d);

for i = N1:N2
    j = i - N1 + 1;
    EjB = conv(E(j, 1:i), Btil);
    
    % Fill G matrix (future control actions)
    if i <= Nu
        G(j, 1:i) = EjB(i:-1:1);
    else
        G(j, :) = EjB(i:-1:i-Nu+1);
    end
    
    % Fill H matrix (past control actions)
    if length(EjB) > i
        H(j, 1:min(nb+d, length(EjB)-i)) = EjB(i+1:end);
    end
end

% GPC gain calculation
Qe = delta * eye(N);
Qu = lambda * eye(Nu);
K_gpc = (G' * Qe * G + Qu) \ (G' * Qe);
K_gpc1 = K_gpc(1, :);  % First row for current control

fprintf('Matrix Dimensions:\n');
fprintf('G: %dx%d, K_gpc: %dx%d\n\n', size(G,1), size(G,2), size(K_gpc,1), size(K_gpc,2));

%% Step 6: Simulation Setup
t = 0:Ts:Tsim;

% nit = length(t);        % Number of iterations
% nin = max(na, nb) + 5;  % Initial samples for settling

nin = 5;
nit = 100 + nin; % número de iterações da simulação

fprintf('=== SIMULATION SETUP ===\n');
fprintf('nit=%d, nin=%d\n\n', nit, nin);

% Initialize arrays with proper dimensions
saidas = zeros(nit, 1);     % System outputs
entradas = zeros(nit, 1);   % Control inputs  
du = zeros(nit, 1);         % Control increments
refs = zeros(nit, 1);       % References

% PWM signals
pwm_a = zeros(nit, 1);
pwm_b = zeros(nit, 1); 
pwm_c = zeros(nit, 1);

% Current references in αβ frame
Iref_alpha = zeros(nit, 1);
Iref_beta = zeros(nit, 1);

% Carrier signal for PWM
carrier = sawtooth(2*pi*f_sw*t, 0.5);

%% Step 7: Reference Generation
% Generate sinusoidal current references
I_ref_peak = 5;  % Peak current reference [A]

for k = 1:nit
    if k > nin
        % Step reference for testing
        refs(k) = I_ref_peak;
        
        % For sinusoidal tracking:
        % refs(k) = I_ref_peak * sin(2*pi*f_fund*t(k));
    end
end

% Transform to αβ frame (for current control)
for k = 1:nit
    Iref_alpha(k) = refs(k);  % For simplicity, use α component only
    Iref_beta(k) = 0;         % Zero β component for balanced operation
end

%% Step 8: Main Simulation Loop
fprintf('Starting GPC + PWM simulation...\n');

for k = nin:nit-1
    %% Step 8.1: System Model (Plant Simulation)
    % Simulate the discrete system
    past_outputs = saidas(max(1,k-na):k);
    past_inputs = entradas(max(1,k-nb-d):k-1);
    
    % Pad with zeros if necessary
    if length(past_outputs) < na+1
        past_outputs = [zeros(na+1-length(past_outputs), 1); past_outputs];
    end
    if length(past_inputs) < nb+d+1
        past_inputs = [zeros(nb+d+1-length(past_inputs), 1); past_inputs];
    end
    
    % System output (simplified model)
    saidas(k+1) = -A(2:end) * past_outputs(end-1:-1:end-na)' + ...
                   B * past_inputs(end:-1:end-nb)';
    
    %% Step 8.2: GPC Control Law
    % Reference vector
    R = refs(k) * ones(N, 1);
    
    % Free response calculation
    f = zeros(N, 1);
    current_output = saidas(k);
    
    % Simplified free response (using current output)
    for j = 1:N
        f(j) = current_output;  % Basic prediction - can be enhanced
    end
    
    % Add contribution from past control actions if H is not empty
    if ~isempty(H) && k > nb+d
        past_du = du(k-1:-1:max(1,k-nb-d));
        if length(past_du) < nb+d
            past_du = [past_du; zeros(nb+d-length(past_du), 1)];
        end
        f = f + H * past_du;
    end
    
    % Control increment calculation
    du(k) = K_gpc1 * (R - f);
    entradas(k) = entradas(k-1) + du(k);
    
    %% Step 8.3: Space Vector Modulation
    V_ref = entradas(k);  % Voltage reference from GPC
    
    % Normalize voltage reference
    V_mag = min(abs(V_ref), Vdc/sqrt(3)) * sign(V_ref);
    theta = 0;  % Assume zero angle for simplicity
    
    % SVM calculations
    sector = 1;  % Simplified - always sector 1
    alpha = 0;
    
    T1 = Ts * V_mag * sin(pi/3 - alpha) / Vdc;
    T2 = Ts * V_mag * sin(alpha) / Vdc;
    T0 = Ts - T1 - T2;
    
    % Duty cycles for sector 1
    Ta = (T1 + T2 + T0/2) / Ts;
    Tb = (T2 + T0/2) / Ts;
    Tc = T0/2 / Ts;
    
    % Apply PWM
    pwm_a(k) = double(Ta > (carrier(k) + 1)/2);
    pwm_b(k) = double(Tb > (carrier(k) + 1)/2);
    pwm_c(k) = double(Tc > (carrier(k) + 1)/2);
    
    % Optional: Add measurement noise
    measurement_noise = 0.01 * randn;
    saidas(k+1) = saidas(k+1) + measurement_noise;
end

fprintf('Simulation completed!\n');

%% Step 9: Results Analysis and Plotting
fprintf('\n=== RESULTS ANALYSIS ===\n');

% Calculate performance metrics
steady_start = find(t > 0.02, 1);
steady_state_error = mean(abs(saidas(steady_start:end) - refs(steady_start:end)));
control_effort = mean(abs(du(nin:end)));

fprintf('Steady-state error: %.4f\n', steady_state_error);
fprintf('Average control effort: %.4f\n', control_effort);

% Plot results
figure('Position', [100, 100, 1200, 800]);

% Subplot 1: Output tracking
subplot(3,2,1);
plot(t, saidas, 'b-', 'LineWidth', 2); hold on;
plot(t, refs, 'r--', 'LineWidth', 2);
title('GPC Output Tracking');
xlabel('Time [s]'); ylabel('Output');
legend('System Output', 'Reference', 'Location', 'best');
grid on;

% Subplot 2: Control input
subplot(3,2,2);
plot(t, entradas, 'g-', 'LineWidth', 2);
title('Control Input');
xlabel('Time [s]'); ylabel('Control Signal');
grid on;

% Subplot 3: Control increments
subplot(3,2,3);
plot(t, du, 'm-', 'LineWidth', 2);
title('Control Increments');
xlabel('Time [s]'); ylabel('\Delta u');
grid on;

% Subplot 4: PWM signals (zoomed)
subplot(3,2,4);
t_zoom = t(1000:1100);
plot(t_zoom, pwm_a(1000:1100), 'r', 'LineWidth', 2); hold on;
plot(t_zoom, pwm_b(1000:1100), 'g', 'LineWidth', 2);
plot(t_zoom, pwm_c(1000:1100), 'b', 'LineWidth', 2);
plot(t_zoom, carrier(1000:1100), 'k--', 'LineWidth', 1);
title('PWM Signals (Zoom)');
xlabel('Time [s]'); ylabel('State');
legend('PWM A', 'PWM B', 'PWM C', 'Carrier', 'Location', 'best');
grid on;
ylim([-0.1 1.1]);

% Subplot 5: Tracking error
subplot(3,2,5);
error = saidas - refs;
plot(t, error, 'r-', 'LineWidth', 2);
title('Tracking Error');
xlabel('Time [s]'); ylabel('Error');
grid on;

% Subplot 6: Performance summary
subplot(3,2,6);
metrics = [steady_state_error, control_effort, max(abs(du)), rms(error(steady_start:end))];
metric_names = {'SS Error', 'Control Effort', 'Max \Delta u', 'RMS Error'};
bar(metrics);
set(gca, 'XTickLabel', metric_names);
title('Performance Metrics');
ylabel('Value');
grid on;

sgtitle('GPC Control with PWM Modulation - 3-Phase VSI', 'FontSize', 14, 'FontWeight', 'bold');

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