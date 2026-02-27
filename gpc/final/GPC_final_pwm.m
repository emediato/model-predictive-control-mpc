%%% Exemplo duplo integrador do capítulo GPC 
%%% cálculo do ganho irrestrito
%%% cálculo das funções de transferência equivalentes.

clear all, 
close all, 
clc

%% Sistema MIMO - Conversor Trifásico GPC

% Parâmetros do hardware
Vd = 60;        % Tensão DC-link [V]
Vdc = Vd;       % Tensão DC-link
R = 22.5;       % Resistência [Ω]
L = 0.1;       % Indutância [H]
fs = 10e3;        % 20 kHz 
f_control_gpc = 1e3;        % 1 kHz

Ts_pwm = 1/fs;
Ts_control = 1/f_control_gpc;
Ts = Ts_control;    % Período de amostragem [s]


T = 0.5;            % Tempo total de simulação
t_control = 0:Ts_control:T;
t_pwm = 0:Ts_pwm:T;

%  CÁLCULO DO TEMPO TOTAL - 2 CICLOS COMPLETOS \n\n');

f = 60;              % Frequência fundamental (Hz)
% Parâmetros do sistema (do item 1)
f_pwm = fs;      % 10 kHz - Frequência PWM
f_control = f_control_gpc;   % 1 kHz - Frequência do controle GPC

Ts_pwm = 1/f_pwm;   % Período PWM = 0.1 ms
Ts_control = 1/f_control; % Período controle = 1 ms
ratio = f_pwm / f_control; % Razão = 10

% Para cada atualização do controle, o PWM executa 10 ciclos. 
% O PWM opera mais rápido para garantir boa resolução de tensão/corrente no motor.
% O duty cycle calculado pelo controlador permanece constante durante esses 10 ciclos de PWM.
% O controle opera mais lentamente para calcular a correção GPC.
% DOUBT A razão ratio = 10 garante uma separação adequada entre Controle (cálculo) Atuação (PWM)

fprintf('PARÂMETROS:\n');
fprintf('Frequência PWM: %.1f kHz → Ts_pwm = %.3f ms\n', f_pwm/1000, Ts_pwm*1000);
fprintf('Frequência Controle: %.1f Hz → Ts_control = %.1f ms\n', f_control, Ts_control*1000);

% Modelo
% Modelo Contínuo do Conversor - Eixos dq
I_2 = eye(2);

% Matriz de Transformação de Clark (αβ → dq)
K = 2/3 * [1 -1/2 -1/2; 
           0 sqrt(3)/2 -sqrt(3)/2];

% Matrizes de espaço de estados
F = -R/L * I_2;  % Matriz de estados
G = (Vd/(2*L)) * K; % Matriz de entrada

fprintf('  MATRIZES DO SISTEMA CONTÍNUO  \n');
disp('Matriz F (dinâmica dos estados):');
disp(F);
disp('Matriz G (entradas):');
disp(G);

% Matrizes discretas
A = expm(F*Ts); % Matriz de estados discreta
B = -F^(-1) * (I_2 - A) * G; % Matriz de entrada discreta

fprintf('\nMatriz A discreta (exp(F*Ts)):');
disp(A);
fprintf('Matriz B discreta:');
disp(B);

s = tf('s');

% Para obter G(1,1): Id/Vd (com Iq = 0)
% Sistema : x1 = G(1,1)*y1 + G(1,2)*y2

% Matriz de transferência contínua
G_s = ss(F, G, I_2, 0); % C = I_2, D = 0
G_tf = tf(G_s);

fprintf('\n  FUNÇÕES DE TRANSFERÊNCIA MIMO  \n');
disp('Matriz de Transferência G(s):');
G_tf

% Modelo perturbação da planta
s = tf('s');
Gq_s = ss(1e-10*F, 1e-10*G, I_2, 0);
Gq_tf = tf(Gq_s);

% FUNÇÕES DE TRANSFERÊNCIA DAS PERTURBAÇÕES 
disp('Matriz de Transferência Gq(s):');
Gq_tf 
% 
%% Discretização
% %%% discretização do processo
Gz_aux = c2d(G_tf,Ts,'zoh');
Gzq = c2d(Gq_tf,Ts,'zoh'); % perts de entrada

%%% adicionando os atrasos discretos manualmente
z  = tf('z',Ts);

delays = [z^-1, z^-0, z^-0;  z^-0, z^-0, z^-0];
% Gz =  Gz_aux.*delays

Gz=Gz_aux;

% Gz = [Gz(1,1), Gz(1,2); Gz(2,1), Gz(2,2);]
% delays = [z^-1, z^-3; z^-1, z^-3;];
delays_q = [z^-1, z^-3;    z^-7, z^-3;];


Gz =  Gz.*delays;


Gqz = Gzq.*delays;
% Gqz =  [Gzq(1,1); Gzq(2,1)];
% Gqz = Gqz.*[z^-8;  z^-4];
sys = minreal(ss([Gz,Gqz]))


% MIMO GPC
m = 3; % número de entradas da planta - manipuladas
n = 2; % número de saidas da planta - controladas
mq = 1; % número de perturbações da planta

A = sys.A;
B = sys.B(:,1:m);
Bq = sys.B(:,m+1:end);
C = sys.C;
Cq = sys.D(:,m+1:end);


%% Sintonia do SSMPC

N1 = [2,3]; %% horizonte de predição inicial - incluir atraso - dimensão 1 x n
N2 = [30,20]; %% horizonte predição 
Ny = N2-N1+1;
Nu = [5,5,5]; %% horizonte de controle - dimensão 1 x m 
Nq = 1; %% horizonte de perturbação - dimensão 1 x mq

lambda = [1,1,1]./Nu; % ponderação do incremento de controle
delta = [5,5]./Ny;  % ponderação dos erros futuros

af = 0; % polo do filtro de referência

Umax = [20, 20, 20];
Umin = [-10,-10, -10];
dumax = [5, 5, 5];
dumin = -dumax;

%% obtenção das matrizes do SSMPC

%%% montando matrizes desconsiderando as diferenças de horizontes, assim,
%%% utiliza-se o valor máximo de N2 e Nu.
temp1 = 0;

F = [];
Ii = [];
for i=1:max(N2)
    temp1 = temp1+A^i;
    
    Ftemp = [C*temp1 -C*temp1];
    F = [F;Ftemp];
    Ii = [Ii; eye(n)];
end


Fq = [C*Bq+Cq, -C*Bq];
temp1 = eye(size(A));
for i=1:max(N2)-1
    temp1 = temp1+A^i;
    
    Fqtemp = [C*temp1*Bq+Cq, -C*temp1*Bq];
    Fq = [Fq;Fqtemp];
    
end



%%% G    
G0 = [];
temp1 = eye(size(A));
for i=1:max(N2)
    
    G0 = [G0;C*temp1*B];
    
    temp1 = temp1 + A^i;
end

G = G0;
for i= 1:max(Nu)-1
    temp = [zeros(i*n,m);G0(1:end-i*n,:)];
    G = [G,temp];
end


%%% Gq
G0q = [Cq];
temp1 = eye(size(A));
for i=1:max(N2)-1
    
    G0q = [G0q;C*temp1*Bq+Cq];
    
    temp1 = temp1 + A^i;
end

Gq = [];
for i=0:Nq-1
    temp = [zeros(i*n,mq);G0q(1:end-i*n,:)];
    Gq = [Gq,temp];
end


%%% matrizes de ponderações
Qei = diag(delta);

Qe = [];
for i=1:max(N2)
    Qe = blkdiag(Qe,Qei);
end

Qui = diag(lambda);

Qu = [];
for i=1:max(Nu)
    Qu = blkdiag(Qu,Qui);
end


%%% Selecionando as parcelas das matrizes de interesse de acordo com os
%%% horizontes estabelecidos na sintonia.

indc = [];
for i=1:m
%     indc = [indc, (1:Nu(i))*m-1+(i-1)];
    indc = [indc, i:m:(Nu(i)-1)*m+i];
end

indc = sort(indc);

indcq = [];
for i=1:mq
    indcq = [indcq, i:mq:(Nq(i)-1)*mq+i];
%     indcq = [indcq, (1:Nq(i))*mq-1+(i-1)];
end

indcq = sort(indcq);


indl = [];
for i=1:n
%     indl = [indl, (N1(i):N2(i))*n-1+(i-1)];
    indl = [indl, (N1(i)-1)*n+i:n:(N2(i)-1)*n+i];
end

indl = sort(indl);


F = F(indl,:);
Ii = Ii(indl,:);
Fq = Fq(indl,:);
G = G(indl,indc);

if(max(Nq)>0)
    Gq = Gq(indl,indcq);
end

Qu = Qu(indc,indc);
Qe = Qe(indl,indl);

Kmpc = (G'*Qe*G+Qu)\(G'*Qe);
Kmpc1 = Kmpc(1:m,:);

Hqp = 2*(G'*Qe*G+Qu);
fqp1 = -2*G'*Qe;

%%% matrizes de restrições
LB = repmat(dumin',max(Nu),1);
LB = LB(indc,:);
UB = repmat(dumax',max(Nu),1);
UB = UB(indc,:);

I = eye(m);
T = [];
for i=1:max(Nu)
    T = [T;
          repmat(I,1,i),zeros(m,max(Nu)*m-i*m)];
end

T = T(indc,indc);

Rbar = [T;
        -T];


%% vetores de simulação
nin = 10; % iteração inicial (para inicialização correta)
nit = nin+150; % número de iterações da simulação

refs = zeros(n,nit); % vetor das referências
perts = zeros(mq,nit); % vetor da perturbação

refs(1,nin+10:end) = 20;
refs(2,nin+50:end)=15;

perts(1,nin+100:end) = 0;

%% simulação do SSMPC
estados = zeros(size(A,1),nit);
saidas = zeros(n,nit); % vetor da saída
entradas = zeros(m,nit); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle

rfant = zeros(n,1);

for k=nin+1:nit
    %%% simulação do processo
    estados(:,k+2) = A*estados(:,k-1)+B*entradas(:,k-1); % +Bq*perts(:,k-1);
    aux = (C*estados(:,k));
              
    saidas(:,k) = aux;

    %%% vetor de referências
    rf = af*rfant + (1-af)*refs(:,k);
    rfant = rf;
    R = repmat(rf,max(N2),1);
    R = R(indl,:);
    
    %%% cálculo da resposta livre;
    f = Ii*saidas(:,k)...
        + F*[estados(:,k);
           estados(:,k-1)];
    
    %%% cálculo das matrizes de restrição
    rbarmax = repmat([Umax'-entradas(:,k-1)],max(Nu),1);
    rbarmin = repmat([-Umin'+entradas(:,k-1)],max(Nu),1);
    rbar = [rbarmax(indc,:);
             rbarmin(indc,:)];
    
    fqp = fqp1*(R-f);
    
    %%% cálculo do incremento de controle ótimo
    % duOti= Kmpc1*(R-f);
    duOti = quadprog(Hqp,fqp,Rbar,rbar,[],[],LB,UB);
    du(1:m,k) = duOti(1:m);
    
    %%% cálculo do sinal de controle ótimo
    entradas(:,k) = du(:,k)+entradas(:,k-1);

    if (k == 60) || (k == 120) || (k == 150)
        i_init = 10; % iteração inicial (para inicialização correta)
        i_sim = i_init + 1000; % número de iterações da simulação
        out = aux; %i_alfa, i_beta

        in = entradas(:,k);
        fprintf('in %f', in)

        [v_abc, duty, t_pwm] = applyPWM(in, out, i_init, i_sim, Ts_pwm, Ts_control, t_control, fs, Vdc, Nu, N1, N2, delta, lambda)
            k_inicio = max(1, i_init);
        
        %entradas(:,k) = v_abc.fase;

        fprintf('fase v_abc  %f', v_abc.fase);
        fprintf('linha v_abc  %f', v_abc.linha);
        fprintf('polo v_abc  %f', v_abc.polo);
        fprintf('duty  %f', duty.sinal);
        fprintf('t_pwm  length %f', length(t_pwm));

        % in = cálculo do sinal de controle ótimo
        % i_init = iteração inicial (para inicialização correta)
        % i_sim = número de iterações da simulação
        % out = saida da simulação do processo
        % Ts_pwm = periodo chaveamento
        % Ts_pwm = periodo controlador
        % t_control = vetor de instante de controle
        % fs = frequencia chaveamento
        % Vdc = tensao barramento
        % Nu = horizonte de controle
        % N1 = horizonte de predição inicial - incluir atraso - dimensão 1 x n
        % N2 = horizonte predição FINAL
        % delta = ponderação do incremento de controle
        % lambda = ponderação dos erros futuros
    end

end


%% Geração das figuras
cores = jet(10);
cores = cores(1:end-1,:);
tamlinha = 4;
tamletra = 12;

hf = figure
h=subplot(3,1,1)
plot(saidas(1,nin+1:nit)','LineWidth',tamlinha,'Color',cores(7,:))
hold on
plot(saidas(2,nin+1:nit)','LineWidth',tamlinha,'Color',cores(1,:))
plot(refs(:,nin+1:nit)','-.','LineWidth',tamlinha,'Color',cores(9,:))
% ylim([0 2])
ylabel('Controladas','FontSize', tamletra)
hl0 = legend('y_1','y_2','Refs.','Location','NorthEast')

title(sprintf(['Horizontes: Nu=[%d %d], N1=[%d %d], N2=[%d %d], ', ...
               'delta=[%.4f %.4f], lambda=[%.4f %.4f]',...
               'Umax = [%d %d]', 'Umax = [%d %d]', 'dumax = [%d %d]'], ...
               Nu(1), Nu(2), N1(1), N1(2), N2(1), N2(2), ...
               delta(1), delta(2), lambda(1), lambda(2), ...
               Umax(1), Umax(2), Umin(1), Umin(2), dumax(1), dumax(2)) );


set(h, 'FontSize', tamletra);
grid on

h = subplot(3,1,2)
plot(entradas(1,nin+1:nit)','LineWidth',tamlinha,'Color',cores(5,:))
hold on
plot(entradas(2,nin+1:nit)','LineWidth',tamlinha,'Color',cores(2,:))
% plot(entradas(3,nin+1:nit)','LineWidth',tamlinha,'Color',cores(3,:))

ylabel('Manipuladas','FontSize', tamletra)
hl1 = legend('u_1','u_2','Location','NorthEast')

grid on
set(h, 'FontSize', tamletra);
% ylim([-5 5])

h = subplot(3,1,3)
plot(du(1,nin+1:nit)','LineWidth',tamlinha,'Color',cores(5,:))
hold on
plot(du(2,nin+1:nit)','LineWidth',tamlinha,'Color',cores(2,:))
% plot(du(3,nin+1:nit)','LineWidth',tamlinha,'Color',cores(3,:))

hl2 = legend('\Delta u_1','\Delta u_2','Location','NorthEast')


ylabel('\Delta u','FontSize', tamletra)
xlabel('Tempo (segundos)','FontSize', tamletra)
grid on
% ylim([-5 5])
set(h, 'FontSize', tamletra);

% h=subplot(3,1,3)
% plot(t,erro(1:nit-N2-1),'LineWidth',tamlinha,'Color',cores(1,:))
% title('Erro', 'FontSize', tamletra);
% xlabel('tempo', 'FontSize', tamletra);
% set(h, 'FontSize', tamletra);


% hl0.Position = [0.8220 0.6325 0.1357 0.1619];
% hl1.Position = [0.8274 0.4933 0.1161 0.1202];
% print('exemploColunaMentanol2','-depsc')


%% Função para após calcular sinal de controle, 
%% aplicar PWM por alguns ciclos
function [v_abc, duty, t_pwm] = applyPWM(in, out, i_init, i_sim, ...
                                           Ts_pwm, Ts_control, t_control, ...
                                           fsw, Vdc, Nu, N1, N2, delta, lambda)

        % applyPWM(in, out, i_init, i_sim, Ts_pwm, Ts_control, t_control, fs, Vdc, Nu, N1, N2, delta, lambda)
        % in = cálculo do sinal de controle ótimo
        % i_init = iteração inicial (para inicialização correta)
        % i_sim = número de iterações da simulação
        % out = saida da simulação do processo
        % Ts_pwm = periodo chaveamento
        % Ts_pwm = periodo controlador
        % t_control = vetor de instante de controle
        % fs = frequencia chaveamento
        % Vdc = tensao barramento
        % Nu = horizonte de controle
        % N1 = horizonte de predição inicial - incluir atraso - dimensão 1 x n
        % N2 = horizonte predição FINAL
        % delta = ponderação do incremento de controle
        % lambda = ponderação dos erros futuros


    % Park (abc → dq) - para medição de correntes
    K_park = @(theta) (2/3) * [cos(theta),   cos(theta - 2*pi/3),   cos(theta - 4*pi/3);
                               -sin(theta), -sin(theta - 2*pi/3), -sin(theta - 4*pi/3)];

    % Clark Transformation Matrix
    K = 2/3 * [1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2];
    % inverse Clarke (αβ -> abc)
    K_inv = [1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2]; % Inverse Clark

    i_alfa = out(1); i_beta = out(2);
    ik = K * [1, 0 ; 0, 1; -1, -1] * [i_alfa; i_beta];
    ia = ik(1);
    ib = ik(2);
    ic= ia-ib;;

    % corrente de fase
    i_abc = [ia ; ib ; ic ];
    fprintf('i_abc %f', i_abc)

    % pré alocação
    t = ((i_init:i_init)-i_init)*Ts_pwm;
    t_pwm = linspace(0, 10*Ts_pwm, length(t));

    ind = i_init:i_sim;

    k_inicio = max(1, i_init);
    k_fim    = min(i_sim, length(t_control));
    
    t_ctrl_janela = t_control(k_inicio:k_fim);   % [1 x n_janela]
    n_janela      = length(t_ctrl_janela);
    
    % Razão de períodos: quantas amostras PWM por ciclo de controle
    % Deve ser inteiro — sincronismo DSP/modulador
    ratio = round(Ts_control / Ts_pwm);
    
    % EIXO DE TEMPO PWM
    % -------------------------------------------------------------------------
    % Resolução fina: passo Ts_pwm dentro da janela de controle.
    % Cada intervalo [k, k+1] do controlador contém exatamente
    % "ratio" amostras PWM — correspondendo a ratio ciclos de portadora.
    
    t_pwm = t_ctrl_janela(1) : Ts_pwm : t_ctrl_janela(end);
    n_pwm = length(t_pwm);
    
    Tsim = 0.02;         % tempo de simulação (s)
    
    fs = 20000;         % taxa de amostragem para simulação dos sinais
    t = 0:1/fs:Tsim;
    f0 = 60; omega = 2*pi*f0;


    % Vetor | Estado     | v_ag        | v_bg         | v_cg
    % ------|------------|-------------|--------------|-------------
    %  V0   | (0, 0, 0)  |      0      |      0       |      0       ← nulo
    %  V1   | (1, 0, 0)  |  +2Vdc/3    |  -Vdc/3      |  -Vdc/3
    %  V2   | (1, 1, 0)  |  +Vdc/3     |  +Vdc/3      |  -2Vdc/3
    %  V3   | (0, 1, 0)  |  -Vdc/3     |  +2Vdc/3     |  -Vdc/3
    %  V4   | (0, 1, 1)  |  -2Vdc/3    |  +Vdc/3      |  +Vdc/3
    %  V5   | (0, 0, 1)  |  -Vdc/3     |  -Vdc/3      |  +2Vdc/3
    %  V6   | (1, 0, 1)  |  +Vdc/3     |  -2Vdc/3     |  +Vdc/3
    %  V7   | (1, 1, 1)  |      0      |      0       |      0       ← nulo

    % carrier triangular (centred -1..1)
    % carrier = sawtooth(2*pi*fsw*t, 0.5); % triangular entre -1 e 1
    carrier = 0.5 + 0.5*sawtooth(2*pi*fsw*t, 0.5);  % Triangular [0,1]
    
%     u_d_array = out(1)*ones(1,(length((ind))))
%     u_q_array = out(2)*ones(1,(length((ind))))

    % u_x = switch_state( find(g==min(g)) , 1:end )' v_alpha =
    % v_abc = 2/3 * (Vdc/2 * in(1) * exp(j*0) + Vdc/2 * in(2) * exp(j*2/3) + Vdc/2 * in(3) * exp(j*4*pi/3) );
    % pré alocação
    ma = zeros(size(t)); mb=ma; mc=ma;
    va_ref = in(1) * ones(1,(length(t))); 
    vb_ref = in(2) * ones(1,(length(t))); 
    vc_ref = in(3) * ones(1,(length(t))); 
    va = in(1); vb = in(2) ; vc = in(3);

    v_alphabeta_ref = K * [in(1); in(2); in(3)];
    
    fprintf('Vd: %.1f a %.1f V\n', min(in), max(in));
    fprintf('Vq: %.1f a %.1f V\n', min(in), max(in));
    fprintf('Vabc %.1f %.1f %.1f \n', va, vb, vc);

    ga = zeros(size(t)); gb=ga; gc=ga;
    idx = 1;
    
    for k=2:ind
        % u_d = u_d_array(k-1);
        % u_q = u_q_array(k-1);
        u_d = v_alphabeta_ref(1);
        u_q = v_alphabeta_ref(2);

        t_start = (k-1)*Ts_control;
        t_end = min(k*Ts_control, Tsim);

        idxs = find(t>=t_start & t < t_end);
    
        for jj = 1:length(idxs)
            ti = t(idxs(jj));
            theta = omega*ti;                 % ângulo instantâneo (ex.: grid angle)
            fprintf('theta %.1f \n', theta);

            % inverse Park
            v_alpha =  cos(theta)*u_d - sin(theta)*u_q;
            v_beta  =  sin(theta)*u_d + cos(theta)*u_q;
    
            % inverse Clarke => va,vb,vc (fase-neutral references)
            vabc = K_inv * [v_alpha; v_beta];
            va = vabc(1); vb = vabc(2); vc = vabc(3);

            % modulation indices  
            ma_i = 0.5*va / Vdc;
            mb_i = 0.5*vb / Vdc;
            mc_i = 0.5*vc / Vdc;
    
            % limit
            ma_i = max(min(ma_i,1),0); duty_a = ma_i;
            mb_i = max(min(mb_i,1),0); duty_b = mb_i;
            mc_i = max(min(mc_i,1),0); duty_c = mc_i;

            % TENSÕES DE POLO (Phase-to-DC-Negative)
            v_AN = duty_a * Vdc;    % [1 x n_pwm], valores em {0, Vdc}
            v_BN = duty_b * Vdc;
            v_CN = duty_c * Vdc;

            % Nolimite daregião linear o modulador sintetiza uma tensão máxima igual a
            % Vao_fase = 2/pi * Vdc;
            % TENSÕES DE FASE (Phase-to-Neutral)
            v_ag = (2*v_AN - v_BN - v_CN) / 3;    % [1 x n_pwm]
            v_bg = (2*v_BN - v_AN - v_CN) / 3;
            v_cg = (2*v_CN - v_AN - v_BN) / 3;

            % TENSÕES DE LINHA (Line-to-Line)
            v_ab = v_AN - v_BN;    % equivalente a v_ag - v_bg
            v_bc = v_BN - v_CN;
            v_ca = v_CN - v_AN;

            % Tensões geradas pelo inversor
            v_abc.fase  = [v_ag;  v_bg;  v_cg ];   % [3 x n_pwm] tensões de fase
            v_abc.linha = [v_ab;  v_bc;  v_ca ];   % [3 x n_pwm] tensões de linha
            v_abc.polo  = [v_AN;  v_BN;  v_CN ];   % [3 x n_pwm] tensões de polo

            % Sinais de chaveamento
            duty.sinal  = [ma_i; mb_i; mc_i];  % [3 x n_pwm] binário

            % store
            va_ref(idxs(jj)) = va; vb_ref(idxs(jj)) = vb; vc_ref(idxs(jj)) = vc;
            ma(idxs(jj)) = ma_i; mb(idxs(jj)) = mb_i; mc(idxs(jj)) = mc_i;
    
            % compare with carrier to make gate signals
            carrier_val = carrier(idxs(jj));
            ga(idxs(jj)) = ma_i > carrier_val;
            gb(idxs(jj)) = mb_i > carrier_val;
            gc(idxs(jj)) = mc_i > carrier_val;

            %va = vabc(1); vb = vabc(2); vc = vabc(3);

            
        end
    end
    
    % plots rápidos
    figure;
    n_plots = 5;
    subplot(n_plots,1,1); 
    plot(t,va_ref,'b', t, vb_ref,'r', t, vc_ref,'g', 'LineWidth',5); 
    ylabel('v^*_ph (V)'); 
    legend('va^*','vb^*','vc^*')
    title(sprintf(['Horizontes: Nu=[%d %d], N1=[%d %d], N2=[%d %d], ', ...
                   'delta=[%.4f %.4f], lambda=[%.4f %.4f]'], ...
                   Nu(1), Nu(2), N1(1), N1(2), N2(1), N2(2), ...
                   delta(1), delta(2), lambda(1), lambda(2)));

    subplot(n_plots,1,2); plot(t,ma,t,mb,t,mc,'LineWidth',5);
    ylabel('mod indices');
    
    subplot(n_plots,1,3); 
    plot(t,ga,'b','LineWidth',1);
    ylabel('gates A (0/1)'); xlabel('time (s)');

    subplot(n_plots,1,4); 
    plot(t,gb,'r','LineWidth',1);
    ylabel('gates B (0/1)'); xlabel('time (s)');

    subplot(n_plots,1,5); 
    plot(t,gc,'g','LineWidth',1);
    %plot(t,ga,'b',t,gb,'r',t,gc,'g','LineWidth',1);
    ylabel('gates C (0/1)'); xlabel('time (s)');
    a = 1;
%    subplot(n_plots,1,4); 
%     plot(t,ga,'b',t,gb,'r',t,gc,'g','LineWidth',1);
%     ylabel('gates (0/1)'); xlabel('time (s)')
end
