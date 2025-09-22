% https://github.com/treezao/Livro_mpc_codes/blob/main/Volume1/Capitulo3-DMC/EstudoDeCaso2/estudo1_dmc.m

% Codigo que implementa DMC
%https://github.com/treezao/Livro_mpc_codes/blob/main/Volume1/Capitulo3-DMC/Exemplo3.22/caso_dmc.m
clc;
clear all;
close all;

% definicao do processo discreto
% y(k)=a1y(k-1)+a2y(k-2)+......a_(na-1)y(k-(na-1))
% +b0u(k-1)+b1u(k-2)+...+b_(nb-1)u(k-(nb-1))

% ex 3.1 pag 49
numG = 20;          % Numerador de G
denG = [1 2 1];    % Denominador de G [s^2, s, termo constante]
G_tf = tf(numG, denG) + 70.04; 

Nss = 70;

t_final = 14; % Choose appropriate final time
t = linspace(0, t_final, Nss); % Nss equally spaced time points




figure(1)
y = step(G_tf, t);           % Step response with 60 points
plot(t, y)                   % Plot the response
grid on
title('Step Response')
xlabel('Time (seconds)')
ylabel('Amplitude')


% Verify we have exactly 60 elements
fprintf('Number of elements in response: %d\n', length(y));
fprintf('Number of time points: %d\n', length(t));


num_coeffs = G_tf.Numerator;
den_coeffs = G_tf.Denominator;

% caldeira temperatura
begin_y = 70.04;
last_y = 90;

%%
rho = 1;
amplitude_step = 10;
delta = amplitude_step;
% 
% qtt = 70;
% g0 = linspace(begin_y, last_y, qtt);
%g0=[0.1000,    0.3000,    0.7000, 0.9];
%g0 = [g0,1,1,1,1,1];

init_valor = 70; %degrees
% u(k) = alfa1 * u (k-1) + alfa2 + ... + beta * y(k) + gamma * r(k);

coefs_g = (y - init_valor ) * 1/delta 
g=coefs_g';

figure;
plot(g);
title('Resposta ao degrau')


%% sistema com atraso
% coeficientes do modelo com atraso
% Add 6 zeros to the start of the vector 'g'
g0_padded_start = [zeros(1, 6), g];

% % Display the first 15 elements to show the change
% disp('Original vector start:');
% disp(g(1:15))
% disp('Padded vector start (6 zeros added):');
% disp(g0_padded_start(1:21))

figure;
plot(g0_padded_start)
title('Resposta ao degrau trocador de calor com atraso')

%WRONG
% coefs_g2 = (g0_padded_start - init_valor ) * 1/delta 
%g2=coefs_g;


%% MPC DMC
% parâmetros de ajuste
N1 = 1; %horizonte de predição inicial
N2 = 20; % horizonte de predição final
N = N2-N1+1; % horizonte de predição


Nu = 5; % horizonte de controle
Nss=40; % horizonte de modelo



% ex 3.3 pagina 56
% horizonte controle
Nu = 5 ;
% horizonte predicao
Nx1 = 1;
Nx2 = 20;

Nx=Nx2;

rbar=1e-1;
% ubar=Nu*rbar % + M*dbar ;
% xbar=Nx*rbar % + Mx*dbar ;

% G tem dimensao 5x20
a = g';

G=zeros(5,5);
% i = 1
G(1:5,1) = a(1:Nu,1);
% i = 2
G(1:5,2) = [zeros(1,1)'; a(1:4,1)];
% i = 3
G(1:5,3) = [zeros(1,2)'; a(1:3,1)];
% i = 4
G(1:5,4) = [zeros(1,3)'; a(1:2,1)];


for i = 1:Nu % (N-1) termos subsequentesi=
    if i==1
        G(1:Nu,i) = a(1:Nu,1);
    else
        aux = Nu-i+1;
        G(1:Nu,i) = [zeros(1,i-1)'; a(1:aux,1)];
    end
end

disp(G);


h = zeros(Nss, Nx);

for i = 1:Nx  
    for j = 1:Nss
        if j+ i < Nss
            h(j,i) = g(j+i) - g(i) ;
        else
             h(j,i) = g(Nss) - g(i) ;
        end
    end
end


% G_correct = zeros(Nx,Nu);
% i=1
% G_correct(1:Nx,i) = a(1:Nx,1)
% i=2
% G_correct(1:Nx,i) = [zeros(1,2)'; a(1:Nx-i+1,1)]
% i=3
% G_correct(1:Nx,i) = [zeros(1,2)'; a(1:Nx-i+1,1)]
% i=4
% G_correct(1:Nx,i) = [zeros(1,3)'; a(1:Nx-i+1,1)]

for i = 1:Nu % (N-1) termos subsequentes
    if i==1
        G_correct(1:Nx,i) = a(1:Nx,1);
    else
        aux = Nx-i+1;
        G_correct(1:Nx,i) = [zeros(1,i-1)'; a(1:aux,1)];
    end
end

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle

%G = G(N1:end,:);
% 
% Qy = delta*eye(N2-N1+1);
% Qu = lambda*eye(Nu);
% 
% Kdmc = inv(G'*Qy*G+Qu)*G'*Qy
% 
% Kdmc1 = Kdmc(1,:);

%%% calcula a matriz H para o cálculo da resposta livre no caso DMC
H1 = [];
H2 = [];


for i=N1(1):N2(1)
    H1 = [H1;coefs_g(i+1:i+Nss)'];
    H2 = [H2;coefs_g(1:Nss)'];
    
end
H = H1-H2;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                REALIZA A SIMULACAO DO PROCESSO DURANTE N AMOSTRAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% inicialização vetores

nin = Nss+1;
% nit = 40 + nin; % número de iterações da simulação

Ts= 1;
nit = round(400/Ts) + nin; % número de iterações da simulação

entradas = 0*ones(nit,1); % vetor o sinal de controle;
du = zeros(nit,1); % vetor de incrementos de controle

saidas = ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,2); % vetor com as perturbações do sistema
perts(nin+round(200/Ts):end) = -.03;

refs = ones(nit,1); % vetor de referências

Cobar = 1.5;

refs(nin+round(50/Ts):end) = 0.1+Cobar;

erro = zeros(nit,1); % vetor de erros

%% simulação com referencia futura
for i = nin:nit
    %% simulação do modelo processo
    saidas(i) = simModelo(saidas(i-1),entradas(i-1)+perts(i-1),Ts);
    
    erro(i) = refs(i)-saidas(i);
    %% Controlador
    
    %%% resposta livre
    f = H*du(i-1:-1:i-Nss) + saidas(i);
    
    %%% referências
    R = ones(N2-N1+1,1)*refs(i);
    
    %%% calculo do incremento de controle ótimo    
    % com restrições
    fqp = fqp1*(R-f);
    
    rbar = [repelem(umax-entradas(i-1),Nu)';
             repelem(entradas(i-1)-umin,Nu)';
             ymax-f;
             -ymin+f];
    [X,FVAL,EXITFLAG] = quadprog(Hqp,fqp,Rbar,rbar,[],[],LB,UB);
    
    %%% caso dê infactível, remover a restrição na saída
    if(EXITFLAG==-2)
        rbar = [repelem(umax-entradas(i-1),Nu)';
                 repelem(entradas(i-1)-umin,Nu)'];
        X = quadprog(Hqp,fqp,Rbar(1:2*Nu,:),rbar,[],[],LB,UB);
    end

    du(i) = X(1);
    
    
    % sem restrições
%     du(i) = Kdmc1*(R-f);
    
    %%% aplicar no sistema
    entradas(i) = entradas(i-1)+du(i);    
end

%% plots
t = ((nin:nit)-nin)*Ts;
vx = nin:nit;

cores = gray(4);
cores = cores(1:end-1,:);


hf = figure
h=subplot(3,1,1)
plot(t,saidas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,refs(vx),'--','LineWidth',tamlinha,'Color',cores(2,:))
ylim([1.3 1.45])
hl = legend('DMC','Referência','Location','NorthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('C_o (mol/L)','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(3,1,2)
plot(t,entradas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
h.YTick = [3.2400 3.2800 3.3200]
ylabel('C_f (mol/L)','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

h = subplot(3,1,3)
plot(t,du(vx),'LineWidth',tamlinha,'Color',cores(1,:))
ylabel('\Delta u','FontSize', tamletra)
xlabel('Tempo (s)','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

hf.Position = tamfigura;
hl.Position = [0.7370 0.6532 0.2054 0.1231];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = 10
len = N;

delta_u = repmat([0]',len-1,1);
h_impulse_delta = repmat([0]',len-1,1);
y_new = repmat([0]',len-1,1);
eta_y = repmat([0]',len-1,1);


for k = 1:N
    %Atualiza o vetor de incremento de controle aplicados no processo com o
    %valor atual calculado

    delta_u(k) = y(k+1) / g(k);

    % resposta ao impulso
    % h(k) = g(k) - g(k-1);
    %h(k) = g(k) - g(k-1)
    h_impulse_delta(k+1) = g(k+1) - g(k); % delta_g

    % resposta do sistema - pagina 47
    y_new(k) = (h_impulse_delta(k) * delta_u(k));

    y(k+1) = y_new(k);

    eta_y(k) = y(k+1) - y(1);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% page 53
% delta_u (k+j-i) = y(k+j) / gi
deulta_u_0 = y(2)/g(1);
deulta_u_1 = y(3)/g(2);


% ^y (k+j|k) = delta_u (k+j-i) * gi + 
% y_chapeu_0 = 




% function [u_fut]=controle_u(Ka,w_ref,fu,D,p,fd)
% %                       Algoritmo do Controle Preditivo DMC
% %                       -----------------------------------
% % [u_fut]=controle_u(Ka,w_ref,fu,D,p,fd)  onde
% %       u_fut    = vetor de controle futuros
% %       Ka       = vetor contendo a primeira linha da matriz K (vem da funçao "matriz_K.m" )
% %       w_ref    = vetor da trajetoria de referencia (vem da funçao "processo.m")
% %       fu       = resposta livre devido ao controle (vem da funçao "resposta_livre.m")
% %       D        = matriz originada da resposta ao salto da perturbacao mensuravel com o 
% %                  numero de linhas igual ao valor do horizonte de prediçao e o numero de 
% %                  colunas igual ao valor do horizonte de controle (vem da funçao "matriz_G.m")
% %       p        = vetor contendo a perturbacao (vem da funçao "processo.m")
% %       fd       = resposta livre devido a perturbacoes mensuraveis (vem da funçao "resposta_livre.m")
% 
% 
% u_fut=Ka*(w_ref-fu-D*p-fd); % Com modelo de perturbação
% % u_fut=Ka*(w_ref-fu);       % Sem modelo de pertrubação
% 
% 
% function [f]=resposta_livre(hor_pred,delta_u,y_medido,g,a_rsalto)
% %                       Algoritmo do Controle Preditivo DMC
% %                       -----------------------------------
% % [f]=resposta_livre(hor_pred,delta_u,y_medido,g,a_rsalto)  onde
% %       f        = vetor com a respota livre do sistema no instante t
% %       hor_pred = horizonte de predição (entrar com número inteiro)
% %       delta_u  = vetor dos valores da variação dos controles passados (delta u(t+i)) (vem da função "processo.m")
% %       y_medido = valor da saída medida no instante t (vem da função "processo.m")
% %       g        = vetor contendo a resposta ao salto da funçao de transferência "g" (vem da função "processo.m")
% %       a_rsalto = Número de amostras que o processo leva para chegar em R.P. (vem da funçao "processo.m")
% %
% %   A cada instante t da simulação do processo esta função é chamada.
% %   O "for" abaixo gera o vetor da resposta livre a cada instante t com
% %   o número de posições igual ao horizonte de predição.
% %   O "for" interno é utilizado para gerar o vetor gk, onde SOMA(g(k+i)-g(i))
% %   com a soma indo de i até N.
% p=hor_pred;  %p ->horizonte de predição
% 
% %DEFINIÇAO DA RESPOSTA LIVRE f
% 
% u_pas1=delta_u(1:a_rsalto,1);
% 
% for t=1:p
% 
%     for e=1:a_rsalto
% 
%        gk(1,e)=(g(t+e)-g(e));
% 
%     end
% 
%     f(t,1)=y_medido+(gk*u_pas1);
% 
% end 

