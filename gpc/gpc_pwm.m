%% CONTROLE GPC + MODULAÇÃO PWM - INVERSOR TRIFÁSICO
% Integração do controle preditivo generalizado com modulação PWM
clear; clc; close all;

%% Parâmetros do Sistema
fprintf('=== CONTROLE GPC + PWM TRIFÁSICO ===\n\n');

% Parâmetros do Inversor
Vdc = 400;              % Tensão DC-link [V]
f_chaveamento = 10e3;   % Frequência de chaveamento [Hz]  
f_fundamental = 60;     % Frequência fundamental [Hz]
Ts = 1/f_chaveamento;   % Período de amostragem [s]
Tsim = 0.1;             % Tempo de simulação [s]

% Parâmetros de Carga
R = 10;                 % Resistência [Ω]
L = 10e-3;              % Indutância [H]
fprintf('Parâmetros do Sistema:\n');
fprintf('  Vdc = %.0f V\n', Vdc);
fprintf('  f_sw = %.1f kHz\n', f_chaveamento/1e3);
fprintf('  R = %.1f Ω, L = %.1f mH\n\n', R, L*1e3);

%% Transformada de Clark
K_clark = 2/3 * [1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2];
K_clark_inv = [1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2]; % Inversa aproximada

%% Modelo de Estado no Referencial αβ (para GPC)
fprintf('Calculando modelo de estado GPC...\n');

I_2 = eye(2);
F = -R/L * I_2;
G = (Vdc/(2*L)) * K_clark(:,1:2); % Considerando apenas 2 fases para controle

% Matrizes de estado discretas
A_d = expm(F*Ts); % matrix exponential  
B_d = -F^(-1) * (I_2 - A_d) * G;

fprintf('Matriz A (discreta):\n'); disp(A_d);
fprintf('Matriz B (discreta):\n'); disp(B_d);

%% Parâmetros GPC
lambda = 0.01;          % Peso do esforço de controle
lambda_i = 10;          % Peso do erro de corrente
N1 = 1;                 % Horizonte inferior
N2 = 5;                 % Horizonte superior  
Nu = 3;                 % Horizonte de controle

fprintf('\nParâmetros GPC:\n');
fprintf('  N1 = %d, N2 = %d, Nu = %d\n', N1, N2, Nu);
fprintf('  lambda = %.3f, lambda_i = %.1f\n', lambda, lambda_i);

%% Cálculo da Matriz G para GPC
fprintf('Calculando matriz G do GPC...\n');

% Dimensões
n_states = size(A_d, 1);
n_inputs = size(B_d, 2);

% Matriz G (resposta ao degrau)
G_gpc = zeros(N2*n_states, Nu*n_inputs);

% Preencher matriz G
for j = 1:Nu
    for i = j:N2
        if i == j
            G_gpc((i-1)*n_states+1:i*n_states, (j-1)*n_inputs+1:j*n_inputs) = B_d;
        else
            temp = eye(n_states);
            for k = j:i-1
                temp = A_d^(i-k) * temp;
            end
            G_gpc((i-1)*n_states+1:i*n_states, (j-1)*n_inputs+1:j*n_inputs) = temp * B_d;
        end
    end
end

% Matriz de ponderação
Q = lambda_i * eye(N2*n_states);  % Ponderação do erro
R = lambda * eye(Nu*n_inputs);    % Ponderação do controle

% Ganho do GPC
K_gpc = (G_gpc' * Q * G_gpc + R) \ (G_gpc' * Q);
K_gpc = K_gpc(1:n_inputs, :); % Apenas primeira linha (controle atual)

fprintf('Dimensão da matriz G: %dx%d\n', size(G_gpc,1), size(G_gpc,2));
fprintf('Dimensão do ganho K_gpc: %dx%d\n', size(K_gpc,1), size(K_gpc,2));

%% Referências de Corrente
fprintf('Gerando referências de corrente...\n');

t = 0:Ts:Tsim;
N = length(t);

% Referência de corrente trifásica desejada
I_ref_peak = 5; % [A] - pico da corrente de referência

Iref_a = I_ref_peak * sin(2*pi*f_fundamental*t);
Iref_b = I_ref_peak * sin(2*pi*f_fundamental*t - 2*pi/3);
Iref_c = I_ref_peak * sin(2*pi*f_fundamental*t + 2*pi/3);

% Transformar para referencial αβ
Iref_alpha = zeros(1, N);
Iref_beta = zeros(1, N);

for k = 1:N
    I_abc = [Iref_a(k); Iref_b(k); Iref_c(k)];
    I_alphabeta = K_clark * I_abc;
    Iref_alpha(k) = I_alphabeta(1);
    Iref_beta(k) = I_alphabeta(2);
end

% Vetor de referência para GPC
I_ref_gpc = zeros(N2*n_states, N);
for k = 1:N
    I_ref_gpc(:, k) = repmat([Iref_alpha(k); Iref_beta(k)], N2, 1);
end

%% Simulação do Sistema com Controle GPC
fprintf('Iniciando simulação GPC + PWM...\n');

% Variáveis de estado
I_alpha = zeros(1, N);
I_beta = zeros(1, N);
I_alpha_meas = zeros(1, N);
I_beta_meas = zeros(1, N);

% Variáveis de controle
V_alpha_ref = zeros(1, N);
V_beta_ref = zeros(1, N);
V_abc_ref = zeros(3, N);
duty_cycles = zeros(3, N);

% Sinais PWM
PWM_A = zeros(1, N);
PWM_B = zeros(1, N);
PWM_C = zeros(1, N);

% Portadora triangular
carrier = sawtooth(2*pi*f_chaveamento*t, 0.5);

% Ruído de medição
noise_gain = 0.1; % Ruído na medição de corrente

% Loop principal de controle
for k = 2:N-1
    %% MEDIÇÃO (com ruído)
    I_alpha_meas(k) = I_alpha(k) + noise_gain * randn;
    I_beta_meas(k) = I_beta(k) + noise_gain * randn;
    
    %% CONTROLE GPC
    if k > N2
        % Estado atual
        current_state = [I_alpha_meas(k); I_beta_meas(k)];
        
        % Vetor de referência futuro
        ref_vector = I_ref_gpc(:, k);
        
        % Predição de estado livre (sem ações futuras)
        free_response = zeros(N2*n_states, 1);
        x_pred = current_state;
        for j = 1:N2
            free_response((j-1)*n_states+1:j*n_states) = x_pred;
            x_pred = A_d * x_pred;
        end
        
        % Erro de predição
        error_pred = ref_vector - free_response;
        
        % Lei de controle GPC
        delta_u = K_gpc * error_pred;
        
        % Tensão de referência em αβ
        V_alpha_ref(k) = delta_u(1);
        V_beta_ref(k) = delta_u(2);
        
        % Saturação
        V_max = Vdc/sqrt(3); % Limite de tensão
        V_alpha_ref(k) = max(min(V_alpha_ref(k), V_max), -V_max);
        V_beta_ref(k) = max(min(V_beta_ref(k), V_max), -V_max);
    end
    
    %% MODULAÇÃO PWM - Transformação para ABC
    if k > 1
        % Transformar tensão αβ para ABC
        V_alphabeta = [V_alpha_ref(k); V_beta_ref(k)];
        V_abc = K_clark_inv * V_alphabeta;
        V_abc_ref(:, k) = V_abc;
        
        % Duty cycles normalizados
        duty_A = 0.5 + V_abc(1) / Vdc;
        duty_B = 0.5 + V_abc(2) / Vdc;
        duty_C = 0.5 + V_abc(3) / Vdc;
        
        % Saturação dos duty cycles
        duty_A = max(min(duty_A, 1), 0);
        duty_B = max(min(duty_B, 1), 0);
        duty_C = max(min(duty_C, 1), 0);
        
        duty_cycles(:, k) = [duty_A; duty_B; duty_C];
        
        % Geração dos sinais PWM
        PWM_A(k) = duty_A > (carrier(k) + 1)/2;
        PWM_B(k) = duty_B > (carrier(k) + 1)/2;
        PWM_C(k) = duty_C > (carrier(k) + 1)/2;
    end
    
    %% MODELO DA PLANTA (Atualização dos estados)
    % Tensões de fase aplicadas
    V_an = (PWM_A(k) - 1/3*(PWM_A(k) + PWM_B(k) + PWM_C(k))) * Vdc;
    V_bn = (PWM_B(k) - 1/3*(PWM_A(k) + PWM_B(k) + PWM_C(k))) * Vdc;
    V_cn = (PWM_C(k) - 1/3*(PWM_A(k) + PWM_B(k) + PWM_C(k))) * Vdc;
    
    % Transformar para αβ
    V_abc_actual = [V_an; V_bn; V_cn];
    V_alphabeta_actual = K_clark * V_abc_actual;
    
    % Atualização das correntes (modelo discreto)
    I_alpha(k+1) = A_d(1,1)*I_alpha(k) + A_d(1,2)*I_beta(k) + ...
                   B_d(1,1)*V_alphabeta_actual(1) + B_d(1,2)*V_alphabeta_actual(2);
    I_beta(k+1) = A_d(2,1)*I_alpha(k) + A_d(2,2)*I_beta(k) + ...
                  B_d(2,1)*V_alphabeta_actual(1) + B_d(2,2)*V_alphabeta_actual(2);
end

fprintf('Simulação concluída!\n');

%% Simulação do Sistema com Controle GPC
fprintf('Iniciando simulação GPC + PWM...\n');

% Variáveis de estado
I_alpha = zeros(1, N);
I_beta = zeros(1, N);
I_alpha_meas = zeros(1, N);
I_beta_meas = zeros(1, N);

% Variáveis de controle
V_alpha_ref = zeros(1, N);
V_beta_ref = zeros(1, N);
V_abc_ref = zeros(3, N);
duty_cycles = zeros(3, N);

% Sinais PWM
PWM_A = zeros(1, N);
PWM_B = zeros(1, N);
PWM_C = zeros(1, N);

% Portadora triangular
carrier = sawtooth(2*pi*f_chaveamento*t, 0.5);

% Ruído de medição
noise_gain = 0.1; % Ruído na medição de corrente

% Loop principal de controle
for k = 2:N-1
    %% MEDIÇÃO (com ruído)
    I_alpha_meas(k) = I_alpha(k) + noise_gain * randn;
    I_beta_meas(k) = I_beta(k) + noise_gain * randn;
    
    %% CONTROLE GPC
    if k > N2
        % Estado atual
        current_state = [I_alpha_meas(k); I_beta_meas(k)];
        
        % Vetor de referência futuro
        ref_vector = I_ref_gpc(:, k);
        
        % Predição de estado livre (sem ações futuras)
        free_response = zeros(N2*n_states, 1);
        x_pred = current_state;
        for j = 1:N2
            free_response((j-1)*n_states+1:j*n_states) = x_pred;
            x_pred = A_d * x_pred;
        end
        
        % Erro de predição
        error_pred = ref_vector - free_response;
        
        % Lei de controle GPC
        delta_u = K_gpc * error_pred;
        
        % Tensão de referência em αβ
        V_alpha_ref(k) = delta_u(1);
        V_beta_ref(k) = delta_u(2);
        
        % Saturação
        V_max = Vdc/sqrt(3); % Limite de tensão
        V_alpha_ref(k) = max(min(V_alpha_ref(k), V_max), -V_max);
        V_beta_ref(k) = max(min(V_beta_ref(k), V_max), -V_max);
    end
    
    %% MODULAÇÃO PWM - Transformação para ABC
    if k > 1
        % Transformar tensão αβ para ABC
        V_alphabeta = [V_alpha_ref(k); V_beta_ref(k)];
        V_abc = K_clark_inv * V_alphabeta;
        V_abc_ref(:, k) = V_abc;
        
        % Duty cycles normalizados
        duty_A = 0.5 + V_abc(1) / Vdc;
        duty_B = 0.5 + V_abc(2) / Vdc;
        duty_C = 0.5 + V_abc(3) / Vdc;
        
        % Saturação dos duty cycles
        duty_A = max(min(duty_A, 1), 0);
        duty_B = max(min(duty_B, 1), 0);
        duty_C = max(min(duty_C, 1), 0);
        
        duty_cycles(:, k) = [duty_A; duty_B; duty_C];
        
        % Geração dos sinais PWM
        PWM_A(k) = duty_A > (carrier(k) + 1)/2;
        PWM_B(k) = duty_B > (carrier(k) + 1)/2;
        PWM_C(k) = duty_C > (carrier(k) + 1)/2;
    end
    
    %% MODELO DA PLANTA (Atualização dos estados)
    % Tensões de fase aplicadas
    V_an = (PWM_A(k) - 1/3*(PWM_A(k) + PWM_B(k) + PWM_C(k))) * Vdc;
    V_bn = (PWM_B(k) - 1/3*(PWM_A(k) + PWM_B(k) + PWM_C(k))) * Vdc;
    V_cn = (PWM_C(k) - 1/3*(PWM_A(k) + PWM_B(k) + PWM_C(k))) * Vdc;
    
    % Transformar para αβ
    V_abc_actual = [V_an; V_bn; V_cn];
    V_alphabeta_actual = K_clark * V_abc_actual;
    
    % Atualização das correntes (modelo discreto)
    I_alpha(k+1) = A_d(1,1)*I_alpha(k) + A_d(1,2)*I_beta(k) + ...
                   B_d(1,1)*V_alphabeta_actual(1) + B_d(1,2)*V_alphabeta_actual(2);
    I_beta(k+1) = A_d(2,1)*I_alpha(k) + A_d(2,2)*I_beta(k) + ...
                  B_d(2,1)*V_alphabeta_actual(1) + B_d(2,2)*V_alphabeta_actual(2);
end

fprintf('Simulação concluída!\n');


%% Transformada inversa para obter correntes ABC
fprintf('Calculando correntes de fase...\n');

I_a = zeros(1, N);
I_b = zeros(1, N);
I_c = zeros(1, N);

for k = 1:N
    I_alphabeta = [I_alpha(k); I_beta(k)];
    I_abc = K_clark_inv * I_alphabeta;
    I_a(k) = I_abc(1);
    I_b(k) = I_abc(2);
    I_c(k) = I_abc(3);
end

% Correntes medidas (com ruído)
I_a_meas = zeros(1, N);
I_b_meas = zeros(1, N);
I_c_meas = zeros(1, N);

for k = 1:N
    I_alphabeta_meas = [I_alpha_meas(k); I_beta_meas(k)];
    I_abc_meas = K_clark_inv * I_alphabeta_meas;
    I_a_meas(k) = I_abc_meas(1);
    I_b_meas(k) = I_abc_meas(2);
    I_c_meas(k) = I_abc_meas(3);
end
