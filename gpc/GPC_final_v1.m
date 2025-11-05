%% LC Filter State-Space Model for GPC
% This script implements the continuous and discrete-time models
% for an LC filter with state-space representation
% transformada clark
K = 2/3 * [1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2] ;

I_2 = eye(2);

%**************************
%                       INIT SYSTEM PARAMETERS
%**************************

% hardware lab measurements
Ts = 1e-2; T_s=Ts;  % Sampling time
     % DC link voltage

i_L_k = 5;       % 5 A
v_alpha_k = 230; % 230 V
v_o_k = 220;     % 220 V
i_g_k = 4.5;     % 4.5 A


Vd = 100; VDC = Vd;Vdc = Vd; V_dc = VDC; % voltage
R = 3.5;  %3.5 % ohms
L = 0.0025 ; % henry

% Example parameters for a typical LC filter
Rf = 0.1;        % Resistance [Ω]
Lf = 1e-3;       % Inductance [mH]
Cf = 10e-6;      % Capacitance [μF]
Ts_filter = 100e-6;   %  [100 μs = 10 kHz]

%% Continuous-Time State-Space Matrices
% state matrix
% Kouros 2009 Modeling
F = - R/L * I_2;
G = (Vd/(2*L)) * K;

Ac = exp(F*Ts);
Bc = -F^(-1)*(I_2 - Ac)*G;

[A,B] = c2dm(Ac,Bc,[],[],Ts,'zoh')

% 
A1 = [-Rf/Lf,  -1/Lf;
       1/Cf,    0   ]
B1 = [Vdc/Lf;  0];
C1 = [1 0];


D1 = 0;
[Bc, Ac] = ss2tf(A1, B1, C1, D1);
Tz = c2d(tf(Bc, Ac), Ts) % Discretized
A = cell2mat(Tz.den); %   A(z-¹)
B = cell2mat(Tz.num); %   B(z-¹)

C = [1 0];


D = 0;
Bq = [0.001]; % numerador do modelo da perturbação do processo

% % state-space system
% sysc = ss(A, B, C, D);
% 
% sysd = ss(A, B, C, D, Ts);
% A = [sysd.A]


%[Ad,Bd] = c2dm(sysc.A,sysc.B,sysc.C,[0 0],Ts,'zoh')


d=0; % atraso em relação à manipulada
dq = 0; % atraso em relação à perturbação

nb = size(B,2)-1; % ordem do polinômio B
na = size(A,2)-1; % ordem do polinômio A
nbq = size(Bq,2)-1; % ordem do polinômio Bq


%% parâmetros de ajuste

N1 = 1; %horizonte de predição inicial
N2 = 10; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 3; % horizonte de controle
Nq = 0; % horizonte da ação antecipativa

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle


%% montagem das matrizes
Btil = conv(B,[zeros(1,d) 1]); % incorporação do atraso no numerador B
Bqtil = conv(Bq,[zeros(1,dq) 1]); % incorporação do atraso no numerador Bq


G = zeros(N,N); % matriz dinâmica G
H = zeros(N,nb+d); % matriz dos incrementos passados de controle 
Gq = zeros(N,N-1); % matris dos incrementos da perturbação futura
Hq = zeros(N,nbq+dq); % matriz dos incrementos passados da perturbação 

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

G_gpc = G;
K_gpc = Kmpc;

fprintf('Dimensão da matriz G: %dx%d\n', size(G_gpc,1), size(G_gpc,2));
fprintf('Dimensão do ganho K_gpc: %dx%d\n', size(K_gpc,1), size(K_gpc,2));

%% **********
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

%% Geração dos Sinais de Referência e Portadora
Tsim = 0.1;             % Tempo de simulação [s]

% inicialização vetores
nin = 5;
nit = 100 + nin; % número de iterações da simulação

% Initialize current prediction variables
% I_alpha(k) = (1 - Rf*Ts/Lf) * I_alpha(k-1) + (V_alpha_ref(k) - Rf*I_alpha(k-1)) * Ts/Lf;
% I_beta(k) = (1 - Rf*Ts/Lf) * I_beta(k-1) + (V_beta_ref(k) - Rf*I_beta(k-1)) * Ts/Lf;

t = 0:Ts:Tsim;
N = length(t);

% Índice de modulação e frequência
m = 0.9;                % Índice de modulação (0 < m <= 1)
f_fundamental = 60;     % Frequência fundamental [Hz]
f_ref = f_fundamental;  % Frequência da referência
f_chaveamento = 10e3;

% Sinais de referência trifásicos
Vref_A = m * sin(2*pi*f_ref*t);
Vref_B = m * sin(2*pi*f_ref*t - 2*pi/3);
Vref_C = m * sin(2*pi*f_ref*t + 2*pi/3);

V_alpha_ref = entradas(k) * cos(2*pi*f_ref*k*Ts);
V_beta_ref = entradas(k) * sin(2*pi*f_ref*k*Ts);

V_ref = V_alpha_ref + 1j*V_beta_ref;

% Portadora triangular
fprintf('Gerando portadora triangular e inicializando vetores...\n');
carrier = sawtooth(2*pi*f_chaveamento*t, 0.5); % Triangular centrada em 0

% Portadoras deslocadas para modulação PWM
carrier_upper = carrier;           % Portadora superior
carrier_lower = -carrier;          % Portadora inferior

% Portadora triangular
carrier = sawtooth(2*pi*f_chaveamento*t, 0.5);

% Ruído de medição
noise_gain = 0.1; % Ruído na medição de corrente

entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle


saidas = zeros(nit, 1);     % System outputs
entradas = zeros(nit, 1);   % Control inputs  
du = zeros(nit, 1);         % Control increments
refs = zeros(nit, 1);       % References


perts = zeros(nit,1); % vetor com as perturbações do sistema
perts(nin+50:end) = 0.5;

refs = 0*ones(nit,1); % vetor de referências
refs(nin+10:end) = 1;

pwm_a = zeros(nit,1);
pwm_b = zeros(nit,1);
pwm_c = zeros(nit,1);


% Current references in αβ frame
Iref_alpha = zeros(nit, 1);
Iref_beta = zeros(nit, 1);

I_ref_peak = 3.3;  % Peak current reference [A]


for k = 1:nit
    if k > nin
        % Step reference for testing
        refs(k) = I_ref_peak;
        
        % For sinusoidal tracking:
        %refs(k) = I_ref_peak * sin(2*pi*f_fund*t(k));
    end
end

% Transform to αβ frame (for current control)
for k = 1:nit
    Iref_alpha(k) = refs(k);  % For simplicity, use α component only
    Iref_beta(k) = 0;         % Zero β component for balanced operation
end

% simulação sem filtro de referência
for k = nin:nit
    % modelo processo, não mexer
    saidas(k) = -A(2:end)*saidas(k-1:-1:k-na) ...
                  +B*entradas(k-d-1:-1:k-1-nb-d) ...
                  +Bq*perts(k-dq:-1:k-nbq-dq);    
    current_output = saidas(k);

    % -- Controlador GPC 
    %%% referencias
    R = refs(k)*ones(N,1);

    % Predição de estado livre (sem ações futuras)
    f = zeros(N, 1);
    % Simplified free response (using current output)
    for j = 1:N
        f(j) = current_output;  % Basic prediction - can be enhanced
    end
    
    %%% cálculo da resposta livre;
    % f = F*saidas(k:-1:k-na);
    
    % if(~isempty(H))
    %     f = f + H*du(k-1:-1:k-nb-d); % parcela dos incrementos de controle
    % end
    %%% cálculo da resposta livre
    %  Calcular f corretamente
    f = zeros(N,1);
    
    %  saídas passadas
    for i = 1:N
        for i = 1:N-1
            f(i) = f(i) + 1 * du(k-1:-1:k-nb-d);
        end
    end
    
    % Contribuição dos incrementos de controle passados
    if ~isempty(H) && (k > nb+d)
        for i = 1:N
            %f(i) = f(i) + H(i,:) * du(k-1:-1:k-nb-d);
            f(i) = f(i) + 1 * du(k-1:-1:k-nb-d);
            
        end
    end

    %% Resolve o problema de otimização
    du(k) = Kmpc1*(R-f);
    
    entradas(k) = entradas(k-1)+du(k);
    
    %% --- Space Vector Modulation (SVM) ---
    V_ref = entradas(k); % desired voltage vector in αβ frame

    % Magnitude and angle
    V_mag = abs(V_ref);
    theta = angle(V_ref); % radians
    theta = mod(theta, 2*pi); % normalize angle

    % Step 2: Identify sector (1 to 6)
    sector = floor(theta / (pi/3)) + 1;

    % Step 3: Calculate T1, T2, T0
    alpha = mod(theta, pi/3);

    T1 = T_s * V_mag * sin(pi/3 - alpha) / V_dc;
    T2 = T_s * V_mag * sin(alpha) / V_dc;
    T0 = T_s - T1 - T2;

    % Step 4: Duty cycles for each phase
    switch sector
        case 1
            Ta = (T1 + T2 + T0/2)/T_s;
            Tb = (T2 + T0/2)/T_s;
            Tc = T0/2/T_s;
        case 2
            Ta = (T1 + T0/2)/T_s;
            Tb = (T1 + T2 + T0/2)/T_s;
            Tc = T0/2/T_s;
        case 3
            Ta = T0/2/T_s;
            Tb = (T1 + T2 + T0/2)/T_s;
            Tc = (T2 + T0/2)/T_s;
        case 4
            Ta = T0/2/T_s;
            Tb = (T1 + T0/2)/T_s;
            Tc = (T1 + T2 + T0/2)/T_s;
        case 5
            Ta = (T2 + T0/2)/T_s;
            Tb = T0/2/T_s;
            Tc = (T1 + T2 + T0/2)/T_s;
        case 6
            Ta = (T1 + T2 + T0/2)/T_s;
            Tb = T0/2/T_s;
            Tc = (T1 + T0/2)/T_s;
    end

    % Store PWM signals
    pwm_a(k) = Ta;
    pwm_b(k) = Tb;
    pwm_c(k) = Tc;

end



%  Current Prediction (Equation 9)
% PREDICT_CURRENT Predicts inductor current at k+1
% Equation (9):
% i_Lαβ(k+1) = (1 - Rf*Ts/Lf)i_Lαβ(k) - Ts/Lf(v_αβ(k) - v_oαβ(k))
% Inputs:
%   i_L_k - Current at time k [A]
%   v_alpha_k - Input voltage αβ at time k [V]
%   v_o_k - Output voltage at time k [V]
%   Rf, Lf - Filter parameters
%   Ts - Sampling time [s]
%
% Output:
%   i_L_next - Predicted current at k+1 [A]

a11 = 1 - Rf*Ts/Lf;
b_coeff = Ts/Lf;

i_L_next = a11 * i_L_k - b_coeff * (v_alpha_k - v_o_k);


% end
% 
% 
% % Predict current
% i_L_next = predict_current(i_L_k, v_alpha_k, v_o_k, Rf, Lf, Ts);
% fprintf('i_L(k) = %.2f A  →  i_L(k+1) = %.4f A\n', i_L_k, i_L_next);
% 
% % Predict voltage
% v_o_next = predict_voltage(v_o_k, i_L_next, i_g_k, Cf, Ts);
% fprintf('v_o(k) = %.2f V  →  v_o(k+1) = %.4f V\n', v_o_k, v_o_next);
% 
% % Using complete state prediction
% x_k = [i_L_k; v_o_k];
% u_k = [v_alpha_k; i_g_k];
% [x_next] = predict_state(x_k, u_k, Rf, Lf, Cf, Ts);
% fprintf('\nUsing state-space model:\n');
% fprintf('x(k+1) = [%.4f; %.4f]\n', x_next(1), x_next(2));
% 
% 
% %% Voltage Prediction (Equation 10)
% function v_o_next = predict_voltage(v_o_k, i_L_next, i_g_k, Cf, Ts)
%     % PREDICT_VOLTAGE Predicts output voltage at k+1
%     %
%     % Equation (10):
%     % v_oαβ(k+1) = v_oαβ(k) + Ts/Cf*i_Lαβ(k+1) - Ts/Cf*i_gαβ(k)
%     %
%     % Inputs:
%     %   v_o_k - Output voltage at time k [V]
%     %   i_L_next - Inductor current at k+1 [A]
%     %   i_g_k - Grid current at time k [A]
%     %   Cf - Filter capacitance [F]
%     %   Ts - Sampling time [s]
%     %
%     % Output:
%     %   v_o_next - Predicted voltage at k+1 [V]
%     
%     b_coeff = Ts/Cf;
%     
%     v_o_next = v_o_k + b_coeff * i_L_next - b_coeff * i_g_k;
% end
% 
% %% Function 5: Complete State Prediction (Equation 11)
%     % PREDICT_STATE Complete state prediction using discrete model
%     %
%     % Equation (11):
%     % [i_Lαβ(k+1)]   [a11  a12] [i_Lαβ(k) ]   [b11  b12] [v_αβ(k) ]
%     % [v_oαβ(k+1)] = [a21  a22] [v_oαβ(k) ] + [b21  b22] [i_gαβ(k)]
%     %
%     % Inputs:
%     %   x_k - State vector at k: [i_Lαβ(k); v_oαβ(k)]
%     %   u_k - Input vector at k: [v_αβ(k); i_gαβ(k)]
%     %   Rf, Lf, Cf - Filter parameters
%     %   Ts - Sampling time [s]
%     %
%     % Outputs:
%     %   x_next - State vector at k+1
%     %   Ad, Bd - Discrete matrices (for reference)
% %     
% %     % Get discrete matrices
% %     [Ad, Bd] = lc_filter_discrete(Rf, Lf, Cf, Ts, 'euler');
% %     
% %     % State prediction
% %     x_next = Ad * x_k + Bd * u_k;
% 
% 
% % GPC Dynamic Matrix
% % fprintf('\n=== GPC Dynamic Matrix ===\n');
% Nu = 2;  % Control horizon
% G = build_gpc_dynamic_matrix(Rf, Lf, Cf, Ts, N, Nu);
% fprintf('G matrix dimensions: %dx%d\n', size(G,1), size(G,2));
% fprintf('G matrix:\n');
% disp(G);


%% plots
t = ((nin:nit)-nin)*Ts;
vx = nin:nit;

cores = gray(3);
cores = cores(1:end-1,:);

tamletra = 12;
tamlinha = 4;
tamfigura=55;
hf = figure
h=subplot(3,1,1)
plot(t,saidas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,refs(vx),'--','LineWidth',tamlinha,'Color',cores(2,:))
ylim([0 1.5])
h.YTick = [0 0.5 1 1.5];
hl = legend('Caso 1','Referência','Location','NorthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(3,1,2)
plot(t,entradas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
h.YTick = [-2 -1 0 1 2]
ylim([-2.5 2])

ylabel('Manipulada','FontSize', tamletra)
grid on

set(h, 'FontSize', tamletra);

h = subplot(3,1,3)
plot(t,du(vx),'LineWidth',tamlinha,'Color',cores(1,:))
% h.YTick = [-2 -1 0 1 2]
% ylim([-2.5 2])
grid on
ylabel('\Delta u','FontSize', tamletra)
xlabel('Tempo (amostras)','FontSize', tamletra)

set(h, 'FontSize', tamletra);

hf.Position = tamfigura;
hl.Position = [0.6952 0.6683 0.2054 0.1242];


%% Função para o cálculo da equação diofantina

function [E,F] = diofantina(A,N1,N2)

    % Cálculo dos polinômios E(z) e F(z)
    
    nA = size(A,2); % ordem de A(z)
    
    % Cálculo de F(z)
    f(1,:) = [1 zeros(1,nA-2)]; % inicializa F
    
    for j=1:N2
        for i=1:nA-2
            f(j+1,i) = f(j,i+1)-f(j,1)*A(i+1);
        end
        f(j+1,nA-1) = -f(j,1)*A(nA);
    end
    
    F = f(1+N1:1+N2,:);
    
    % Cálculo de E(z)
    E = zeros(N2);
    e(1) = 1;
    E(1,1) = e(1);
    
    for i=2:N2
        e(i) = f(i,1);
        E(i,1:i) = e;
    end
    
    E = E(N1:N2,:);
end
