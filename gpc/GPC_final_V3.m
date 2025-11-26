%% ========================================================================
%  IMPLEMENTAÇÃO GPC COM PWM PARA SISTEMA TRIFÁSICO
%  Sistema: Conversor CC-CC controlando correntes em coordenadas dq
%% ========================================================================

clear all; close all; clc;

%% 1. PARÂMETROS DO SISTEMA FÍSICO
%==========================================================================
% Parâmetros do conversor trifásico
Vdc = 400;              % Tensão do barramento CC [V]
L = 5e-3;               % Indutância de fase [H]
R = 0.1;                % Resistência de fase [Ω]
C = 1000e-6;            % Capacitância [F]

% Parâmetros elétricos
f_rede = 60;            % Frequência da rede [Hz]
omega_e = 2*pi*f_rede;  % Velocidade angular elétrica [rad/s]

% Frequências de amostragem e PWM
fs = 10e3;              % Frequência de amostragem/controle [Hz]
Ts = 1/fs;              % Período de amostragem [s]
f_pwm = 10e3;           % Frequência PWM [Hz]
Tpwm = 1/f_pwm;         % Período PWM [s]

% Tempo de simulação
t_sim = 0.5;            % Duração da simulação [s]
N_samples = t_sim/Ts;   % Número de amostras

%% 2. MODELO EM ESPAÇO DE ESTADOS (COORDENADAS dq)
%==========================================================================
% Modelo contínuo em dq:
% d/dt [id]   [[-R/L,  omega_e] [id]     [1/L  0 ] [vd]
%      [iq] =  [-omega_e, -R/L]] [iq]  +  [0  1/L]] [vq]
%
% Estados: x = [id, iq]'
% Entradas: u = [vd, vq]'

% Matrizes do sistema contínuo
Ac = [-R/L,     omega_e;
      -omega_e, -R/L];
  
Bc = [1/L,  0;
      0,    1/L];

Cc = eye(2);  % Medimos id e iq diretamente
Dc = zeros(2,2);

% Sistema contínuo
sys_c = ss(Ac, Bc, Cc, Dc);

% Discretização (ZOH - Zero Order Hold)
sys_d = c2d(sys_c, Ts, 'zoh');
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

fprintf('=== MODELO DISCRETIZADO ===\n');
fprintf('Período de amostragem: %.2f us\n', Ts*1e6);
disp('Matriz Ad:'); disp(Ad);
disp('Matriz Bd:'); disp(Bd);

%% 3. PROJETO DO CONTROLADOR GPC
%==========================================================================
% Parâmetros do GPC
N1 = 1;                 % Horizonte mínimo
N2 = 10;                % Horizonte máximo de predição
Nu = 5;                 % Horizonte de controle

% Pesos do GPC
lambda = 0.1;           % Peso do esforço de controle
Q_gpc = eye(2);         % Peso do erro de rastreamento (id, iq)

% Restrições
u_max = 0.95*Vdc/2;     % Tensão máxima dq (índice de modulação)
u_min = -u_max;
du_max = 50;            % Variação máxima de tensão [V]

%% 4. CONSTRUÇÃO DAS MATRIZES DO GPC
%==========================================================================
% Matriz G (dinâmica do sistema)
G = zeros(N2, Nu);
for j = 1:N2
    for i = 1:min(j, Nu)
        if i == 1
            G(j,i) = Cd * Ad^(j-1) * Bd(1,1);  % Para id (primeira coluna de Bd)
        else
            G(j,i) = Cd * Ad^(j-i) * Bd(1,1);
        end
    end
end

% Expandir para 2 saídas (id e iq) e 2 entradas (vd e vq)
G_full = kron(G, eye(2));  % Kronecker para sistemas MIMO

% Matriz H (Hessiana)
H = 2*(G_full' * kron(eye(N2), Q_gpc) * G_full + lambda * eye(2*Nu));

% Garantir H simétrica e positiva definida
H = (H + H')/2;
H = H + 1e-6*eye(size(H));  % Regularização numérica

%% 5. CONFIGURAÇÃO DAS RESTRIÇÕES
%==========================================================================
% Restrições de amplitude: u_min ≤ u ≤ u_max
% Em forma matricial: A_ineq * ΔU ≤ b_ineq

% Matriz de soma acumulada Su
Su = tril(ones(Nu, Nu));
Su_full = kron(Su, eye(2));

% Restrições de amplitude
A_u = [Su_full; -Su_full];

% Restrições de variação
A_du = [eye(2*Nu); -eye(2*Nu)];

% Matriz A completa
A_ineq = [A_u; A_du];

%% 6. INICIALIZAÇÃO DA SIMULAÇÃO
%==========================================================================
% Estados iniciais
x = [0; 0];  % id(0) = 0, iq(0) = 0

% Referências em dq (constantes para este exemplo)
id_ref = 10;   % Corrente d desejada [A]
iq_ref = 5;    % Corrente q desejada [A]

% Vetores para armazenar resultados
t_vec = (0:N_samples-1)*Ts;
id_vec = zeros(1, N_samples);
iq_vec = zeros(1, N_samples);
vd_vec = zeros(1, N_samples);
vq_vec = zeros(1, N_samples);
u_last = [0; 0];  % Ação de controle anterior

fprintf('\n=== INICIANDO SIMULAÇÃO ===\n');

%% 7. LOOP PRINCIPAL DE CONTROLE
%==========================================================================
for k = 1:N_samples
    
    % Estados atuais
    id_vec(k) = x(1);
    iq_vec(k) = x(2);
    
    % Erro de rastreamento
    e = [id_ref; iq_ref] - x;
    
    % Resposta livre (predição sem ação de controle)
    f = zeros(2*N2, 1);
    x_free = x;
    for j = 1:N2
        x_free = Ad * x_free;
        f((j-1)*2+1:j*2) = Cd * x_free;
    end
    
    % Vetor c (gradiente da função objetivo)
    r_vec = kron(ones(N2,1), [id_ref; iq_ref]);  % Vetor de referências
    c = 2*G_full' * kron(eye(N2), Q_gpc) * (f - r_vec);
    
    % Atualizar vetor b das restrições (dependente do estado atual)
    b_u_upper = kron(ones(Nu,1), [u_max; u_max]) - kron(ones(Nu,1), u_last);
    b_u_lower = kron(ones(Nu,1), u_last) - kron(ones(Nu,1), [u_min; u_min]);
    b_du_upper = kron(ones(Nu,1), [du_max; du_max]);
    b_du_lower = kron(ones(Nu,1), [du_max; du_max]);
    
    b_ineq = [b_u_upper; b_u_lower; b_du_upper; b_du_lower];
    
    % Resolver QP: min 0.5*ΔU'*H*ΔU + c'*ΔU  s.t. A_ineq*ΔU ≤ b_ineq
    options = optimoptions('quadprog', 'Display', 'off');
    [dU_opt, ~, exitflag] = quadprog(H, c, A_ineq, b_ineq, [], [], [], [], [], options);
    
    if exitflag ~= 1
        warning('QP não convergiu no instante k=%d', k);
        dU_opt = zeros(2*Nu, 1);
    end
    
    % Aplicar apenas o primeiro elemento (princípio do horizonte deslizante)
    du = dU_opt(1:2);
    u_atual = u_last + du;
    
    % Saturação adicional (segurança)
    u_atual = max(min(u_atual, [u_max; u_max]), [u_min; u_min]);
    
    % Armazenar ações de controle
    vd_vec(k) = u_atual(1);
    vq_vec(k) = u_atual(2);
    
    % Atualizar sistema (simulação da planta)
    x = Ad*x + Bd*u_atual;
    
    % Atualizar controle anterior
    u_last = u_atual;
    
    % Progresso
    if mod(k, floor(N_samples/10)) == 0
        fprintf('Progresso: %.0f%%\n', 100*k/N_samples);
    end
end

fprintf('=== SIMULAÇÃO CONCLUÍDA ===\n');

%% 8. GERAÇÃO DE SINAIS PWM (PÓS-PROCESSAMENTO)
%==========================================================================
% Índices de modulação (normalizado por Vdc/2)
m_d = vd_vec / (Vdc/2);
m_q = vq_vec / (Vdc/2);

% Conversão dq → αβ → abc para visualização
theta = omega_e * t_vec;  % Ângulo elétrico

% Transformação inversa de Park (dq → αβ)
v_alpha = m_d.*cos(theta) - m_q.*sin(theta);
v_beta = m_d.*sin(theta) + m_q.*cos(theta);

% Transformação inversa de Clarke (αβ → abc)
va = v_alpha;
vb = -0.5*v_alpha + (sqrt(3)/2)*v_beta;
vc = -0.5*v_alpha - (sqrt(3)/2)*v_beta;

%% 9. PLOTAGEM DOS RESULTADOS
%==========================================================================
figure('Name', 'Resultados GPC-PWM', 'Position', [100 100 1200 800]);

% Correntes dq
subplot(3,2,1);
plot(t_vec, id_vec, 'b', 'LineWidth', 1.5); hold on;
plot(t_vec, id_ref*ones(size(t_vec)), 'r--', 'LineWidth', 1);
grid on;
xlabel('Tempo [s]'); ylabel('i_d [A]');
title('Corrente i_d (eixo direto)');
legend('i_d', 'Referência');

subplot(3,2,2);
plot(t_vec, iq_vec, 'b', 'LineWidth', 1.5); hold on;
plot(t_vec, iq_ref*ones(size(t_vec)), 'r--', 'LineWidth', 1);
grid on;
xlabel('Tempo [s]'); ylabel('i_q [A]');
title('Corrente i_q (eixo quadratura)');
legend('i_q', 'Referência');

% Tensões de controle dq
subplot(3,2,3);
plot(t_vec, vd_vec, 'LineWidth', 1.5);
grid on;
xlabel('Tempo [s]'); ylabel('v_d [V]');
title('Tensão de controle v_d');

subplot(3,2,4);
plot(t_vec, vq_vec, 'LineWidth', 1.5);
grid on;
xlabel('Tempo [s]'); ylabel('v_q [V]');
title('Tensão de controle v_q');

% Índices de modulação abc
subplot(3,2,5);
plot(t_vec, va, 'r', t_vec, vb, 'g', t_vec, vc, 'b', 'LineWidth', 1);
grid on;
xlabel('Tempo [s]'); ylabel('Índice de modulação');
title('Sinais de modulação ABC (normalizado)');
legend('Fase A', 'Fase B', 'Fase C');
ylim([-1.1 1.1]);

% Trajetória no plano dq
subplot(3,2,6);
plot(id_vec, iq_vec, 'b', 'LineWidth', 1.5); hold on;
plot(id_ref, iq_ref, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('i_d [A]'); ylabel('i_q [A]');
title('Trajetória no plano dq');
legend('Trajetória', 'Referência');
axis equal;

%% 10. ANÁLISE DE DESEMPENHO
%==========================================================================
fprintf('\n=== ANÁLISE DE DESEMPENHO ===\n');
fprintf('Erro RMS em id: %.4f A\n', rms(id_vec - id_ref));
fprintf('Erro RMS em iq: %.4f A\n', rms(iq_vec - iq_ref));
fprintf('Tensão máxima vd: %.2f V\n', max(abs(vd_vec)));
fprintf('Tensão máxima vq: %.2f V\n', max(abs(vq_vec)));
fprintf('Índice de modulação máximo: %.4f\n', max([abs(va), abs(vb), abs(vc)]));

%% 11. SIMULAÇÃO PWM DETALHADA (OPCIONAL)
%==========================================================================
% Para visualizar o sinal PWM real, precisamos sobreamostrar

% Selecionar um intervalo curto para visualização
t_pwm_start = 0.1;  % [s]
t_pwm_end = 0.11;   % [s]
idx_start = find(t_vec >= t_pwm_start, 1);
idx_end = find(t_vec <= t_pwm_end, 1, 'last');

% Oversampling para PWM
oversample = 100;
t_pwm_detail = linspace(t_vec(idx_start), t_vec(idx_end), ...
                        (idx_end-idx_start+1)*oversample);

% Interpolação dos índices de modulação
ma_detail = interp1(t_vec(idx_start:idx_end), va(idx_start:idx_end), t_pwm_detail, 'linear');

% Geração do carrier (triangular)
carrier = sawtooth(2*pi*f_pwm*t_pwm_detail, 0.5);  % Triangular simétrica

% Comparação PWM
pwm_signal = double(ma_detail > carrier);

% Plot PWM
figure('Name', 'Detalhe PWM', 'Position', [100 100 1000 600]);

subplot(2,1,1);
plot(t_pwm_detail*1000, ma_detail, 'b', 'LineWidth', 2); hold on;
plot(t_pwm_detail*1000, carrier, 'r--', 'LineWidth', 1);
grid on;
xlabel('Tempo [ms]'); ylabel('Amplitude');
title('Modulação PWM - Fase A');
legend('Sinal modulante (m_a)', 'Portadora triangular');
ylim([-1.2 1.2]);

subplot(2,1,2);
stairs(t_pwm_detail*1000, pwm_signal*Vdc - Vdc/2, 'LineWidth', 1.5);
grid on;
xlabel('Tempo [ms]'); ylabel('Tensão [V]');
title('Sinal PWM gerado - Fase A');
ylim([-Vdc/2-50 Vdc/2+50]);

fprintf('\n=== SIMULAÇÃO PWM CONCLUÍDA ===\n');
fprintf('Razão PWM/Controle: %.0f\n', f_pwm/fs);
