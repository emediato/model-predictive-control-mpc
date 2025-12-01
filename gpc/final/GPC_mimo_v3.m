clear; clc; close all;

%%
% Entradas (y): [Vd; Vq] - Tensões de controle nos eixos dq
% Saídas (x): [Id; Iq] - Correntes nos eixos dq
% Perturbações (z): [Var_Tensão_Rede; Var_Carga]

%% Sistema MIMO - Conversor Trifásico GPC
clear; clc; close all;

% Parâmetros do hardware
Vd = 60;        % Tensão DC-link [V]
Vdc = Vd;       % Tensão DC-link
R = 22.5;       % Resistência [Ω]
L = 0.10;       % Indutância [H]
fs = 20e3;        % 20 kHz 
f_control_gpc = 1e3;        % 1 kHz

Ts_pwm = 1/fs;
Ts_control = 1/f_control_gpc;
Ts = Ts_control;    % Período de amostragem [s]


T = 0.5;            % Tempo total de simulação
t_control = 0:Ts_control:T;
t_pwm = 0:Ts_pwm:T;

%  CÁLCULO DO TEMPO TOTAL - 2 CICLOS COMPLETOS \n\n');

% Parâmetros do sistema (do item 1)
f_pwm = fs;      % 10 kHz - Frequência PWM
f_control = f_control_gpc;   % 1 kHz - Frequência do controle GPC

Ts_pwm = 1/f_pwm;   % Período PWM = 0.1 ms
Ts_control = 1/f_control; % Período controle = 1 ms

ratio = f_pwm / f_control; % Razão = 10

fprintf('PARÂMETROS:\n');
fprintf('Frequência PWM: %.1f kHz → Ts_pwm = %.3f ms\n', f_pwm/1000, Ts_pwm*1000);
fprintf('Frequência Controle: %.1f Hz → Ts_control = %.1f ms\n', f_control, Ts_control*1000);


fprintf('1. Cálculo do GPC + Transformações: t_calc\n');
fprintf('2. Aplicação do PWM por %d ciclos: %d × Ts_pwm\n', ratio, ratio);
fprintf('3. Medição + Processamento: t_measure\n\n');

% Tempos típicos (estimativas realistas)
t_calc = 0.2e-3;    % 0.2 ms - cálculo GPC + transformações
t_measure = 0.1e-3; % 0.1 ms - medição + processamento

fprintf('TEMPOS ESTIMADOS:\n');
fprintf('Cálculo GPC + Transformações: %.1f ms\n', t_calc*1000);
fprintf('Aplicação PWM (%d ciclos): %.1f ms\n', ratio, ratio*Ts_pwm*1000);
fprintf('Medição + Processamento: %.1f ms\n', t_measure*1000);

% Tempo de UM ciclo completo
T_um_ciclo = t_calc + (ratio * Ts_pwm) + t_measure;

fprintf('\nTEMPO DE UM CICLO COMPLETO:\n');
fprintf('T_um_ciclo = %.1f + %.1f + %.1f = %.1f ms\n', ...
        t_calc*1000, ratio*Ts_pwm*1000, t_measure*1000, T_um_ciclo*1000);

T_total_2_ciclos = 2 * T_um_ciclo;

fprintf('TEMPO TOTAL (2 ciclos completos):\n');
fprintf('T_total = 2 × %.1f ms = %.1f ms\n', T_um_ciclo*1000, T_total_2_ciclos*1000);
fprintf('T_total = %.3f segundos\n\n', T_total_2_ciclos);


f = 60;              % Frequência fundamental (Hz)
Vdc = 60;            % Tensão DC-link


% Referência de operação
x_ref = [0 5]';  % [Id_ref; Iq_ref]

% Definição das variáveis MIMO:
fprintf('  VARIÁVEIS DO SISTEMA MIMO  \n');
fprintf('ENTRADAS (y - manipuladas):\n');
fprintf('  y1: Tensão no eixo d (Vd)\n');
fprintf('  y2: Tensão no eixo q (Vq)\n\n');

fprintf('SAÍDAS (x - controladas):\n');
fprintf('  x1: Corrente no eixo d (Id)\n');
fprintf('  x2: Corrente no eixo q (Iq)\n\n');

fprintf('PERTURBAÇÕES (z - não manipuladas):\n');
fprintf('  z1: Variação da tensão da rede\n');
fprintf('  z2: Variação de carga\n');

%% Modelo
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
% Sistema SISO: x1 = G(1,1)*y1 + G(1,2)*y2

% Matriz de transferência contínua
G_s = ss(F, G, I_2, 0); % C = I_2, D = 0
G_tf = tf(G_s);

fprintf('\n  FUNÇÕES DE TRANSFERÊNCIA MIMO  \n');
disp('Matriz de Transferência G(s):');
G_tf

% Extrair elementos individuais
G11 = G_tf(1,1); % Id/Vd
G12 = G_tf(1,2); % Id/Vq  
G21 = G_tf(2,1); % Iq/Vd
G22 = G_tf(2,2); % Iq/Vq



%% Modelagem das Perturbações Plausíveis
% Perturbação 1: Variação da tensão da rede (±1%)
% Perturbação 2: Variação de carga (mudança de R)

% Matriz de perturbações Gq
Gq_s = ss(F, [0.001*G(:,1), 0.015*G(:,2)], I_2, 0);
Gq_tf = tf(Gq_s);

% FUNÇÕES DE TRANSFERÊNCIA DAS PERTURBAÇÕES 
disp('Matriz de Transferência Gq(s):');
Gq_tf

G=G_tf;
Gq=Gq_tf;

% %%% discretização do processo
Gz = c2d(G,Ts,'zoh');
Gzq = c2d(Gq,Ts,'zoh'); % perts de entrada
% 

atrasoEntrada = totaldelay(Gz);


% MIMO GPC
m = 2; % número de entradas da planta - manipuladas
n = 2; % número de saidas da planta - controladas
mq = 1; % número de perturbações da planta

s = tf('s');


if(mq>0)
    atrasoPert = totaldelay(Gzq);
end
%% Preparação para GPC MIMO
% Horizonte de predição e controle
Np = 20; % Horizonte de predição
Nc = 5;  % Horizonte de controle

% Matrizes de ponderação
Q = diag([1, 1]); % Ponderação das saídas
R = diag([0.1, 0.1]); % Ponderação das entradas

% Extrair matrizes de espaço de estados discretas
sys_d = ss(Gz, 'min');
[Ad, Bd, Cd, Dd] = ssdata(sys_d);
sys_dq = ss(Gzq, 'min');
[Bdq, Cdq, Ddq] = ssdata(sys_dq); % Matriz de perturbações

fprintf('\nMatriz Ad (discreta):\n');
disp(Ad);
fprintf('Matriz Bd (discreta):\n');
disp(Bd);

fprintf('\n  CONFIGURAÇÃO GPC MIMO  \n');
fprintf('Horizonte de predição (Np): %d\n', Np);
fprintf('Horizonte de controle (Nc): %d\n', Nc);
fprintf('Matriz de ponderação Q:\n');
disp(Q);
fprintf('Matriz de ponderação R:\n');
disp(R);



%% Controlador GPC - MIMO
%% PARAMETROS DE AJUSTE DO CONTROLADOR

Nu = [5 5]; %horizontes de controle - dimensão 1 x m
N1 = [9 10]; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = [10 20]; % fim dos horizontes de predição - dimensão 1 x n
N = N2-N1+1; %horizonte de predição - não editar
Nq = [0]; % horizonte de predição da perturbação - dimensão 1 x mq

delta = [0.5 0.9]./N; %ponderação nos erros - dimensão 1 x n
lambda = [0.5 0.9]./Nu; %ponderação nas ações de controle - dimensão 1 x m       

C{1} = [1 -0.92];%% definição do polinômio C
C{2} = [1 -0.92]; 

% fprintf('   - Controlador roda em frequência mais baixa\n');
% fprintf('   - PWM roda em frequência mais alta\n');
% fprintf('   - Razão típica: 10:1 (PWM:Controlador)\n\n');



%% definição do cenário de simulação

nit = 500; %tempo de simulação
nin = 50; %número da iteração inicial

nit = nit+nin;

ind = nin:nit;
t = ((nin:nit)-nin)*Ts;


fs_control = f_control_gpc;   % Frequência do controle GPC
t=0:1/f_control_gpc:nit;

fs_pwm = f_pwm;      % Frequência do PWM

Am_reference = 5;    % Amplitude de referência
% Geração das referências trifásicas
omega = 2*pi*f;
Va_ref = Am_reference * sin(omega*t_control);
Vb_ref = Am_reference * sin(omega*t_control - 2*pi/3);
Vc_ref = Am_reference * sin(omega*t_control - 4*pi/3);

% Transformação Park das referências
d_ref = zeros(size(t_control));
q_ref = zeros(size(t_control));

for k = 1:length(t_control)
    theta_k = omega * t_control(k);
    
    % Transformação Park direta (abc → dq)
    d_ref(k) = (2/3) * (Va_ref(k)*cos(theta_k) + Vb_ref(k)*cos(theta_k - 2*pi/3) + Vc_ref(k)*cos(theta_k - 4*pi/3));
    q_ref(k) = (2/3) * (-Va_ref(k)*sin(theta_k) - Vb_ref(k)*sin(theta_k - 2*pi/3) - Vc_ref(k)*sin(theta_k - 4*pi/3));
end

dq_ref = [d_ref; q_ref];  % [Id_ref; Iq_ref]


refs = dq_ref(1) * ones(n,nit);
perts = dq_ref(2) * ones(mq,nit+max(Nq));

%%%%%%%%%%%% degrau de referencia

Am_reference = 10;    % Amplitude de referência
% Geração das referências trifásicas
omega = 2*pi*f;
Va_ref = Am_reference * sin(omega*t_control);
Vb_ref = Am_reference * sin(omega*t_control - 2*pi/3);
Vc_ref = Am_reference * sin(omega*t_control - 4*pi/3);

% Transformação Park das referências
d_ref = zeros(size(t_control));
q_ref = zeros(size(t_control));

for k = 1:length(t_control)
    theta_k = omega * t_control(k);
    
    % Transformação Park direta (abc → dq)
    d_ref(k) = (2/3) * (Va_ref(k)*cos(theta_k) + Vb_ref(k)*cos(theta_k - 2*pi/3) + Vc_ref(k)*cos(theta_k - 4*pi/3));
    q_ref(k) = (2/3) * (-Va_ref(k)*sin(theta_k) - Vb_ref(k)*sin(theta_k - 2*pi/3) - Vc_ref(k)*sin(theta_k - 4*pi/3));
end

dq_ref = [d_ref; q_ref];  % [Id_ref; Iq_ref]


refs(100,nin+100:end) = dq_ref(1);
refs(200,nin+800:end) = dq_ref(2);
perts(1,nin+130:end) = 2;

%%
%%% obtenção da representação MFD

[Amfd,Bmfd] = MFDredux(Gz.den,Gz.num);

Bqmfd = Bmfd(:,m+1:end);
Bmfd = Bmfd(:,1:m);

%%% montagem das matrizes de predição
F = {};
H = {};
G = {};

nH = zeros(1,m); % qt de dados passados necessários das entradas
nF = zeros(1,n); % qt de dados passados necessários das saídas.

for i=1:n
    
    [Ei,Fi] = diofantina(conv(Amfd{i},[1 -1]),N1(i),N2(i));
    
    F{i} = Fi;
    nF(i) = size(Fi,2);
    
    for j=1:m
        % incorporação do atraso no numerador B
        Btil = conv(Bmfd{i,j}(2:end),[zeros(1,atrasoEntrada(i,j)) 1]); 
        
        Gtemp = [];
        Htemp = [];
        for k=N1(i):N2(i)
            EjB = conv(Ei(k-N1(i)+1,1:k),Btil);
            
            Gtemp(k-N1(i)+1,1:k) = EjB(k:-1:1);    
            Htemp(k-N1(i)+1,:) = EjB(k+1:end);
        end
        G{i,j} = Gtemp(:,1:Nu(j));
        H{i,j} = Htemp;
        
        nhi = size(Htemp,2);
        if(nhi>nH(j))
            nH(j)=nhi;
        end
    end
        
end

%%% ajuste dos tamanhos das matrizes Hi e Hqi para concatenação

for i=1:n
    for j=1:m
        nHi = size(H{i,j},2);
        if(nHi<nH(j))
            H{i,j} = [H{i,j} zeros(N(i),nH(j)-nHi)];
        end
    end
    
end

H = cell2mat(H);
G = cell2mat(G);

Ftemp = [];
for i=1:n
    Ftemp = blkdiag(Ftemp,F{i});
end
F = Ftemp;
   

        
%%% matrizes de ponderação
Qe = diag(repelem(delta,1,N)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle

% H = 2*(G'*Qe*G+Qu);
% invH = inv(H); 

Kdmc = (G'*Qe*G+Qu)\G'*Qe;
Kdmc1 = [];
for i=1:m
    Kdmc1 = [Kdmc1; Kdmc(sum(Nu(1:i-1))+1,:)];
end

%%% inicialização dos vetores para o cálculo da resposta livre sem correção


%% inicializacao dos estados e variaveis da simulação - GPC
% - obtencao dos numerados e denominadores
numR = cell(n,m);
denR = cell(n,m);
numDR = cell(n,mq);
denDR = cell(n,mq);


for i=1:n
	% obtencao das respostas ao degrau das entradas
	for j=1:m
		% obtencao dos numeradores e denominadores para o simulador
        [numR{i}{j}, denR{i}{j}] = tfdata(Gz(i,j),'v');
	end
    
	% obtencao das respostas ao degrau das perts
    for j=1:mq
        [numDR{i}{j},denDR{i}{j}]=tfdata(Gzq(i,j),'v');
	end
end


estadosEntr = cell(n,m);
for i=1:n
    for j=1:m
		estadosEntr{i}{j} = zeros(1,nit);
    end
end

estadosPert = cell(n,mq);
for i=1:n
	for j=1:mq
		estadosPert{i}{j}= zeros(1,nit);
	end    
end

saidas = zeros(n,nit); % vetor da saída
entradas = zeros(m,nit); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle

for k = nin:nit
    %% -- simulador, não alterar
    %%% Simulador MIMO que utiliza estados
    for i=1:n
        aux = 0;
        for j=1:m
            sizeN = size(numR{i}{j},2);
            sizeD = size(denR{i}{j},2);

            str1 = numR{i}{j}*entradas(j,k-atrasoEntrada(i,j):-1:k+1-sizeN-atrasoEntrada(i,j))';
            str2 = -denR{i}{j}(1,2:end)*estadosEntr{i}{j}(k-1:-1:k-sizeD+1)';
            estadosEntr{i}{j}(k) = str1+str2;

            aux = aux + estadosEntr{i}{j}(k);
        end

        for j=1:mq
            sizeN = size(numDR{i}{j}, 2);
            sizeD = size(denDR{i}{j}, 2);

            str1 = numDR{i}{j}*perts(j,k-atrasoPert(i,j):-1:k+1-sizeN-atrasoPert(i,j))';
            str2 = -denDR{i}{j}(1,2:end)*estadosPert{i}{j}(k-1:-1:k-sizeD+1)';
            estadosPert{i}{j}(k) = str1+str2;

            aux = aux + estadosPert{i}{j}(k);
        end
        saidas(i,k) = aux ;%;+ 0.01*randn(1,1);
    end
    %% -- Controlador GPC 
    %%% referencias
    R = [];
    for i=1:n
        R = [R;
             repelem(refs(i,k),N(i),1)];
    end 
    
    %%% calculo da resposta livre
    du1 = [];
    for i=1:m
        du1 = [du1;du(i,k-1:-1:k-nH(i))'];
    end
    
    y1 = [];
    for i=1:n
        y1 = [y1;saidas(i,k:-1:k-nF(i)+1)'];
    end
    
    f = F*y1 + H*du1;


    %% Resolve o problema de otimização
    du(:,k) = Kdmc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end

saidasSFF = saidas;
entradasSFF = entradas;
duSFF = du;




%% GPC com polinômio C 
%% Controlador GPC - MIMO
%%% montagem das matrizes de predição
F = {};
H = {};
G = {};

nH = zeros(1,m); % qt de dados passados necessários das entradas
nF = zeros(1,n); % qt de dados passados necessários das saídas.
nC = zeros(1,n); % qt de dados passados necessários do poli C

for i=1:n
    
    nC(i) = size(C{i},2);
    
    [Ei,Fi] = diofantinaC(conv(Amfd{i},[1 -1]),C{i},N1(i),N2(i));
    
    F{i} = Fi;
    nF(i) = size(Fi,2);
    
    for j=1:m
        % incorporação do atraso no numerador B
        Btil = conv(Bmfd{i,j}(2:end),[zeros(1,atrasoEntrada(i,j)) 1]); 
        
        Gtemp = [];
        Htemp = [];
        for k=N1(i):N2(i)
            EjB = conv(Ei(k-N1(i)+1,1:k),Btil);

            [Mi,Ni] = diofantinaC(C{i},EjB,k,k);
    
            Gtemp(k-N1(i)+1,1:k) = Mi(k:-1:1);    
            Htemp(k-N1(i)+1,:) = Ni;
        end
        G{i,j} = Gtemp(:,1:Nu(j));
        H{i,j} = Htemp;
        
        nH(i,j) = size(H{i,j},2);
    end
end

G = cell2mat(G);

Kgpc = (G'*Qe*G+Qu)\G'*Qe;
Kgpc1 = [];
for i=1:m
    Kgpc1 = [Kgpc1; Kgpc(sum(Nu(1:i-1))+1,:)];
end



%% configuração da simulação GPC com polinômio C

estadosEntr = cell(n,m);
for i=1:n
    for j=1:m
		estadosEntr{i}{j} = zeros(1,nit);
    end
end

estadosPert = cell(n,mq);
for i=1:n
	for j=1:mq
		estadosPert{i}{j}= zeros(1,nit);
	end    
end

saidas = zeros(n,nit); % vetor da saída
entradas = zeros(m,nit); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle

saidasC = zeros(n,nit);

for i=1:n
    for j=1:m
        duC{i,j} = zeros(1,nit);
    end
end


for k = nin:nit
    %% -- simulador, não alterar
    %%% Simulador MIMO que utiliza estados
    for i=1:n
        aux = 0;
        for j=1:m
            sizeN = size(numR{i}{j},2);
            sizeD = size(denR{i}{j},2);

            str1 = numR{i}{j}*entradas(j,k-atrasoEntrada(i,j):-1:k+1-sizeN-atrasoEntrada(i,j))';
            str2 = -denR{i}{j}(1,2:end)*estadosEntr{i}{j}(k-1:-1:k-sizeD+1)';
            estadosEntr{i}{j}(k) = str1+str2;

            aux = aux + estadosEntr{i}{j}(k);
        end

        for j=1:mq
            sizeN = size(numDR{i}{j}, 2);
            sizeD = size(denDR{i}{j}, 2);

            str1 = numDR{i}{j}*perts(j,k-atrasoPert(i,j):-1:k+1-sizeN-atrasoPert(i,j))';
            str2 = -denDR{i}{j}(1,2:end)*estadosPert{i}{j}(k-1:-1:k-sizeD+1)';
            estadosPert{i}{j}(k) = str1+str2;

            aux = aux + estadosPert{i}{j}(k);
        end
        saidas(i,k) = aux ;%;+ 0.01*randn(1,1);
    end
    %% -- Controlador GPC 
    %%% referencias
    R = [];
    for i=1:n
        R = [R;
             repelem(refs(i,k),N(i),1)];
    end 
    
    %%% calculo da resposta livre
    yc = [];
    f = [];
    for i=1:n
        saidasC(i,k) = saidas(i,k) -C{i}(2:end)*saidasC(i,k-1:-1:k-nC(i)+1)'; 
        
        ftemp = F{i}*saidasC(i,k:-1:k-nF(i)+1)';
        
        for j=1:m
            ftemp = ftemp + H{i,j}*duC{i,j}(1,k-1:-1:k-nH(i,j))';
        end
        f = [f;ftemp];
    end

    
    %% Resolve o problema de otimização
    du(:,k) = Kgpc1*(R-f);
    
    for i=1:n
        for j=1:m
           duC{i,j}(1,k) = -C{i}(2:end)*duC{i,j}(k-1:-1:k-nC(i)+1)' + du(j,k);
        end
    end
    
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end

saidasCFF = saidas;
entradasCFF = entradas;
duCFF = du;

%% após calcular sinal de controle, aplicar PWM por alguns ciclos

t = ((nin:nit)-nin)*Ts;
ind = nin:nit;

% Inicializar sinais modulantes
m_abc = zeros(3, length(t_control)); % [ma; mb; mc]

fprintf('Vd: %.1f a %.1f V\n', min(entradas(1,ind)), max(entradas(1,ind)));
fprintf('Vq: %.1f a %.1f V\n', min(entradas(2,ind)), max(entradas(2,ind)));


Vdc = 60;            % tensão do barramento
fsw = 10000;         % freq. PWM (Hz)
Ts_pwm = 1/fsw;
Tsim = 0.04;         % tempo de simulação (s)

fs = 200000;         % taxa de amostragem para simulação dos sinais

t = 0:1/fs:Tsim;
f0 = 60; omega = 2*pi*f0;

% carrier triangular (centred -1..1)
%carrier = sawtooth(2*pi*fsw*t, 0.5); % triangular entre -1 e 1
carrier = 0.5 + 0.5*sawtooth(2*pi*fsw*t, 0.5);  % Triangular [0,1]

u_d_array = saidas(1,ind);
u_q_array = saidas(2,ind);

% inverse Clarke (αβ -> abc)
K_inv = [1, 0; -1/2, sqrt(3)/2; -1/2, -sqrt(3)/2];

% pré alocação
ma = zeros(size(t)); mb=ma; mc=ma;
va_ref = zeros(size(t)); vb_ref=va_ref; vc_ref=va_ref;
ga = zeros(size(t)); gb=ga; gc=ga;

idx = 1;


for k=1:ind
    u_d = u_d_array(k);
    u_q = u_q_array(k);
    t_start = (k-1)*Ts_control;
    t_end = min(k*Ts_control, Tsim);
    idxs = find(t>=t_start & t < t_end);

    for jj = 1:length(idxs)
        ti = t(idxs(jj));
        theta = omega*ti;                 % ângulo instantâneo (ex.: grid angle)
        % inverse Park
        v_alpha =  cos(theta)*u_d - sin(theta)*u_q;
        v_beta  =  sin(theta)*u_d + cos(theta)*u_q;

        % inverse Clarke => va,vb,vc (fase-neutral references)
        vabc = K_inv * [v_alpha; v_beta];
        va = vabc(1); vb = vabc(2); vc = vabc(3);

        % modulation indices (SPWM)
        ma_i = 0.5*va / Vdc;
        mb_i = 0.5*vb / Vdc;
        mc_i = 0.5*vc / Vdc;

        % limit
        ma_i = max(min(ma_i,1),0);
        mb_i = max(min(mb_i,1),0);
        mc_i = max(min(mc_i,1),0);

        % store
        va_ref(idxs(jj)) = va; vb_ref(idxs(jj)) = vb; vc_ref(idxs(jj)) = vc;
        ma(idxs(jj)) = ma_i; mb(idxs(jj)) = mb_i; mc(idxs(jj)) = mc_i;

        % compare with carrier to make gate signals
        carrier_val = carrier(idxs(jj));
        ga(idxs(jj)) = ma_i > carrier_val;
        gb(idxs(jj)) = mb_i > carrier_val;
        gc(idxs(jj)) = mc_i > carrier_val;
    end
end

% plots rápidos
figure;
subplot(3,1,1); 
plot(t,va_ref,'b', t, vb_ref,'r', t, vc_ref,'g', 'LineWidth',5); 
ylabel('v^*_ph (V)'); 
legend('va^*','vb^*','vc^*')
subplot(3,1,2); plot(t,ma,t,mb,t,mc,'LineWidth',5);
ylabel('mod indices');

subplot(3,1,3); 
plot(t,ga,'b',t,gb,'r',t,gc,'g','LineWidth',1);
ylabel('gates (0/1)'); xlabel('time (s)')
title(sprintf(['Horizontes: Nu=[%d %d], N1=[%d %d], N2=[%d %d], ', ...
               'delta=[%.4f %.4f], lambda=[%.4f %.4f]'], ...
               Nu(1), Nu(2), N1(1), N1(2), N2(1), N2(2), ...
               delta(1), delta(2), lambda(1), lambda(2)));

% 

%%

x=2;
tamlinha=5;
cores = lines(4);

l=length(ind);

hf= figure


h=subplot(x,1,1)
plot(t(1:l),saidasSFF(1,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t(1:l),saidasCFF(1,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
legend('y1','y2'); grid on;

plot(t(1:l),saidasSFF(2,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
plot(t(1:l),saidasCFF(2,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))

cores = summer(4);
plot(t(1:l),refs(1,ind),'-.','Color',[0 0 0],'LineWidth',tamlinha-3,'Color',cores(3,:))

plot(t(1:l),refs(2,ind),'-.','Color',[0 0 0],'LineWidth',tamlinha-3,'Color',cores(4,:))

legend('C_l(z^{-1})=1','C_l(z^{-1})\neq 1', 'Refs',  'Location','SouthEast')

ta1 = annotation('textarrow');
ta1.String = 'y_1';
ta2 = annotation('textarrow');
ta2.String = 'y_2';

title(sprintf(['Horizontes: Nu=[%d %d], N1=[%d %d], N2=[%d %d], ', ...
               'delta=[%.4f %.4f], lambda=[%.4f %.4f]'], ...
               Nu(1), Nu(2), N1(1), N1(2), N2(1), N2(2), ...
               delta(1), delta(2), lambda(1), lambda(2)));


tamletra=12
ylabel('Controladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


h=subplot(x,1,2)
cores = jet(4);
plot(t(1:l),entradasSFF(1,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t(1:l),entradasCFF(1,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))

plot(t(1:l),entradasSFF(2,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
plot(t(1:l),entradasCFF(2,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))

legend('C_l(z^{-1})=1','C_l(z^{-1})\neq 1','Location','SouthEast')

ta3 = annotation('textarrow');
ta3.String = 'u_1';

ta4 = annotation('textarrow');
ta4.String = 'u_2';


ylabel('Manipuladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

xlabel('Tempo (segundos)')
% 
% hf.Position = tamfigura;
% hl.Position = [0.7680 0.4753 0.2179 0.2387];
%y1
ta1.Position = [0.2661 0.7084 0.0629 0.0697];
%y2
ta2.Position = [0.39475 0.7919 0.0471 -0.0516];

%u1
ta3.Position = [0.25382 0.3381 -0.0529 0.0568];
%u2
ta4.Position = [0.3961 0.10084 0.0629 0.0697];


%%
% plots rápidos
figure;
x=3;
tamlinha=2
subplot(x,1,1);
plot(ga, 'LineWidth',tamlinha,'Color',cores(1,:)); hold on;
xlabel('k'), ylabel(' switch A ')

subplot(x,1,2);
plot(gb, 'LineWidth',tamlinha,'Color',cores(2,:)); hold on;
xlabel('k'), ylabel(' switch B ')

subplot(x,1,3);
plot(gc, 'LineWidth',tamlinha,'Color',cores(3,:)); hold on;
xlabel('k'), ylabel(' switch C ')


% print('tanque_gpc_polic_q','-depsc')


function [Ao,Bo] = MFDredux(A,B)
n = size(A,1);
m = size(A,2);
tol = 1e-6; % tolerância para indicar raízes iguais.

Ao = {};
Bo = {};
for i=1:n
    raizes = [];
    for j=1:m
        r = roots(A{i,j});
        while(length(r)>0)
            r0 = r(1); % raíz do polinômio
            mult = 0; % multiplicidade
            I = find(abs(r-r0) < tol); % encontrando os índices de raízes iguais com certa tolerância

            mult = size(I,1); % multiplicidade é o número de raízes encontradas
            
            r(I) = []; % removendo as raízes iguais processadas
            
            % verifica se raíz r0 já existia em outras FTs e se a
            % multiplicidade é maior ou menor
            if(length(raizes)>0)
                I = find(abs(raizes(:,1)-r0) < tol);

                if(size(I,1)==0) % se não existia
                    raizes = [raizes; r0 mult];
                else % se existia
                    if(raizes(I,2) < mult) % verifica multiplicidade
                        raizes(I,2) = mult;
                    end % se multiplicade já era maior, não faz nada
                end
            else
                raizes = [r0 mult];
            end
        end
    end
    
    Ao{i} = poly(repelem(raizes(:,1), raizes(:,2)));    
    
    
    for j=1:m
        [q,r] = deconv(Ao{i},A{i,j});
        
        Bo{i,j} = conv(B{i,j},q);
        
    end
end


end



%% Função para o cálculo da equação diofantina

function [E,F] = diofantinaC(A,C,N1,N2)

% Cálculo dos polinômios E(z) e F(z)
% delta = [1 -1]; % delta = 1-z^{-1}
% AD = conv(A,delta) ; % AD = A(z)*Delta(z)

nA1 = size(A,2); % numero de elemnetos de A(z)
nC = size(C,2); % numero de elementos de C(z)


%%% calculo da ordem máxima do polinômio F
nF1 = max(nC,N1-1+nA1)-N1;
nF2 = max(nC,N2-1+nA1)-N2;
nF = max(nF1,nF2);

% Cálculo de F(z)

if(nC>nA1)
    A = [A zeros(1,nC-nA1)];
    nA = size(A,2);

else
    nA = nA1;
end

f(1,:) = [C zeros(1,nA-nC)]; % inicializa F



for j=1:N2
    for i=1:nA-1
        f(j+1,i) = f(j,i+1)-f(j,1)*A(i+1);
    end
%     f(j+1,nA-1) = -f(j,1)*A(nA);
end

F = f(1+N1:1+N2,1:nF);



% Cálculo de E(z)
E = zeros(N2);
e(1) = C(1);
E(1,1) = e(1);

for i=2:N2
    e(i) = f(i,1);
    E(i,1:i) = e;
end

E = E(N1:N2,:);


 
end

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

