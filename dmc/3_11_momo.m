
s = tf('s');

Ts = 0.05; % Periodo de amostragem do processo em segundos

G = 4/(s+2)^2*exp(-10*Ts*s); % Modelo por função de transferência do processo

Gz = c2d(G,Ts,'zoh'); % modelo discretizado
num = Gz.num{1}; % numerador do modelo discreto
den = Gz.den{1}; % denominador do modelo discreto
na = size(den,2)-1; % ordem do denominador
nb = size(num,2)-2; % ordem do numerador
dd = Gz.inputdelay; % atraso discreto

%% parâmetros de ajuste
% definição dos horizontes

% Sempre N1=1 e N2 é definido como 'p' de predição, tal que p < N 
% p=40;
N1 = 11; %horizonte de predição inicial
N2 = 40; % horizonte de predição final
hor_pred=N2;p=N2;

% horizonte de controle 'Nu' ou 'm', sempre m<=p
% m=5; 
Nu = 5; % horizonte de controle


% N2 < d + Nu
delta = 1; % ponderação do erro futuro
lambda = 1; % λ = ponderação do esforço de controle 

Nss=80; % horizonte de modelo
a_rsalto = Nss;  % número de amostras da resposta ao salto 

Gcoef = step(G,Ts:Ts:Nss*Ts); % coeficientes da resposta ao degrau

figure;
step(G,Ts:Ts:Nss*Ts); 
%% montando as matrizes do DMC recursivo

G = zeros(N2,Nu);
% method 1 for G definition
G(:,1) = Gcoef(1:N2,1);

for i=2:Nu
    G(i:end,i) = G(1:end-(i-1),1);    
end

G = G(N1:end,:);

fprintf('G = ')
disp(G);

a = Gcoef;

%% method 2 for G definition for DMC irrestrito 
% for i = 1:Nu % (N-1) termos subsequentesi=
%     if i==1
%         G(1:Nu,i) = a(1:Nu,1);
%     else
%         aux = Nu-i+1;
%         G(1:Nu,i) = [zeros(1,i-1)'; a(1:aux,1)];
%     end
% end
% disp(G);


%%
% DMC 
% ponderaçoes de controle são constante
% metodo 1
sizeG = size(G); sizeG = sizeG(2);
K = inv(G'*G + lambda*eye(sizeG))*G';
Kc=K(1,:);

% metodo 2
I=eye(size(G'*G)); %Matriz identidade
L=lambda*I;    %Lambida x Matriz Identidade
K=(inv((G'*G)+L))*G';  %Matriz K
Ka=K(1,:);

% DMC RECURSIVO
% metodo 3
% page 78
% ajuste de lambda! + λ - ganho k=j
Qe = delta*eye(N2-N1+1);
Qu = lambda*eye(Nu);

Kdmc = inv(G'*Qe*G + Qu) * G'*Qe ;
Kb = Kdmc(1,:);


%% inicialização vetores

% Se Nss = 80 é o horizonte do modelo em amostras
% Tempo total do horizonte = 80 * 0.05 = 4.0 segundos
% Número de amostras = 80 / 0.05 = 1600 amostras

% pagina 91
% Ts = 3
% Horizonte do modelo: Nss * Ts = 80 × 3s = 240 segundos
% Tempo mínimo de simulação: 3 * 240 = 720 segundos (12 minutos)
% 5 * 240 = 1200
Ts=3;
t_final_fast = 3 * Nss * Ts;  % 720s = 12 min
t_final_complet = 4 * Nss * Ts;  % 960s = 16 min
t_final_extensive = 5 * Nss * Ts;  % 1200s = 20 min

t_sim = t_final_extensive + Nss ; % % número de iterações da simulação

%referencias conhecidas
refs = 0*ones(nit,1); % vetor de referências
refs(nin+50:nit) = 100;
refs(nin+250:nit) = 200;
refs(nin+450:nit) = 300;
refs(nin+650:nit) = 200;
refs(nin+850:end) = 100;
% r = [ r(k+N1), r(k+N1+1), r(k+N1+2), ... , r(k+N2))  ] 

% referencias desconhecida
% r = 1_{N2-N1+1} * r(k) = 1_N*r(k)


entradas = 0*ones(Nss+1,1); % vetor o sinal de controle

y_medido = 0 ; %valor da saída medida no instante t 


%       D        = matriz originada da resposta ao salto da perturbacao mensuravel com o 
%                  numero de linhas igual ao valor do horizonte de prediçao e o numero de 
%                  colunas igual ao valor do horizonte de controle (vem da funçao "matriz_G.m")

perts = zeros(nit,1); % vetor com as perturbações do sistema
perts(nin+150:end) = 0;

refs = 0*ones(t_sim,1); % vetor de referências


% definição do vetor dos incrementos de controle passados para resposta livre
u = zeros(1,t_sim);


delta_u(Nss,1)=0;      %vetor com os valores da variacao dos controle passados

h_impulse_delta = repmat([0]',N-1,1);

eta = zeros(t_sim,1); %η =  vetor de erros de prediçao

y_free  = ones(Nss,1)*0; % 0 é o valor inicial da saída do sistema

p=hor_pred; % vetor contendo a perturbacao p ->horizonte de predição


y_new = repmat([0]', N-1,1);  
y = 0*ones(nit,1); % vetor com as saídas do sistema

input = 0*ones(nit,1); % vetor o sinal de controle

%% dmc 

for k = Nss+1 : t_sim-N2-1
  
    % %Calcula a saida no instante t
    % y_medido=-Atio(1,2:na)*yna(1:na-1,1)+Bu_aux(1,2:nb)*unb(atraso+1:nb+atraso-1,1)...
    %     +Bp_aux(1,2:nb)*pnb(atraso_pert+1:nbp+atraso_pert-1,1);
    % y medido!
    y(k) = -den(2:end) * y(k-1 : -1 : k-na) + num * (input(k-dd:-1 : k-nb-dd-1) + perts(k-dd: -1 :k-dd-nb-1));
    
    
    error = refs(k) - y(k);

    % resposta livre
        % DEFINIÇAO DA RESPOSTA LIVRE f resposta_livre(hor_pred,delta_u,y_medido,g,a_rsalto)
    % Calcula a respota livre devido ao controle e a perturbcao
    %Resposta livre com relação à entrada de controle e com y(k)
    u_pas1 = delta_u(1:Nss,1);
    
    % f -> o número de posições igual ao horizonte de predição.
    for t=1:N2+1
        for e=1:N2 % Número de amostras que o processo leva para chegar em R.P. 
           gk(1,e)=(Gcoef(t+e) - Gcoef(e));
        end
        gk(1,80)=(Gcoef(79) - Gcoef(80));
        
        f(t,1)= y(k) + (gk * u_pas1);
    end 

    fu = f + eta;

     %%% referências
    R = ones(N2-N1+1,1)*refs(k);
    
    % incremento do controle otimo

    % u_fut = vetor de controle futuros
    % u_fut=Ka*(R - fu - D * p - y_free(k)); % Com modelo de perturbação
    u_fut = Ka*(R-fu);       % Sem modelo de pertrubação


    delta_u(k)

    input(i) = uAtual;


end


