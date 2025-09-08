% Codigo que implementa DMC

clc;
clear all;
close all;

% definicao do processo discreto
% y(k)=a1y(k-1)+a2y(k-2)+......a_(na-1)y(k-(na-1))
% +b0u(k-1)+b1u(k-2)+...+b_(nb-1)u(k-(nb-1))

% os 'ai' no vetor A
% os 'bi' no vetor B
% A e B servem para o calculo da simulação do processo

Ts=1;
d=3;
A = [1.6375 -0.6703]; 
B = [zeros(1,d) 0.035 0.0307]; % B tem d colunas com zeros por causa do atraso
num = [0.035 0.0307];

%Parametros para calcular o modelo do processo, pode se colocar erro de
%modelagem definindo o atraso e A e B difernete do modelo
dp=4;
Ap = A;
nump = num;
Bp=[zeros(1,dp) nump]; 

den = [1 -A zeros(1,d-1)]; %  nesta forma de colocar o modelo o atraso equivale a d polos em z=0, isso
% aumenta o num de colunas. Alem disso temos que adicionar o '1' no começo do vetor é para garantir o termo 1*y(k)
% dado que A so tem os coeficientes de z^{-1} em diante

P=tf(num,den,Ts);
zpk(P)

%%


% calcula o numero de coeficientes dos polinomios usados para simulacao

na=size(Ap,2); % 'na' recebe o numero de colunas de Ap (neste exemplo = 2 colunas)
nb=size(Bp,2); % 'nb' recebe o numero de colunas de Bp (neste exemplo = 6 colunas)


% definicao dos vetores onde se guardam as informações passadas da saída y e entradas u

y_ant=zeros(1,na); % numero de colunas de na
u_ant=zeros(1,nb);


% N = número de pontos da resposta ao degrau

N=30;
g=step(P,N);


% definição do vetor dos incrementos de controle passados para resposta livre

u_livre=zeros(1,N);


% definição dos horizontes
% Sempre N1=1 e N2 é definido como 'p' de predição, tal que p <N 
p=20;

% horizonte de controle 'Nu' ou 'm', sempre m<=p
m=5;

% calcula matriz G

G=zeros(p,m); % inicializa a matriz G só com zeros, p linhas e m colunas
G(:,1)=g(1:p); % a primeira coluna da matriz G recebe as amostras da resposta ao degrau
%pegou-se p amostras(20), menos que o vetor N, 30 elementos


for i=2:m % da segunda coluna até a m-ésima culuna
    ga=[zeros(1,i-1) g(1:p-i+1)']; % ga' é o vetor coluna que pega o valor da coluna à esquerda e da um shift com '0' pra baixo
    G(:,i)=ga';
end

% define simulação
pontos=220;

% define parametros do DMC
lambda=10;


% Calculo da Matriz Mn

Mn=(G'*G + lambda*eye(m)); 

% Inverte Mn, multiplica por G' e pega a primeira linha

Mn1=inv(Mn)*G'; % Mn1 é a matriz de ganho k, deltaU = k(w-f)
qn=Mn1(1,:); % pega só a primeira linha de Mn1 ( k1 )


% define vetores de pontos elementos para simular
% yp do proceso, y predição, u controle r=ref 

% controle para graficos
u=zeros(1,pontos);

%vetor de saidas

yp=zeros(1,pontos); %y processo
ypr=zeros(1,pontos); %y processo contabilizando a perturbação

% def referencia

w=[zeros(1,10) 1*ones(1,pontos-10)]; % vetor com 220 pontos, os 10 primeiros são '0' e o restante são '1'

% coloca filtro 1 ordem

af=0. ;% parametro para ajuste do filtro de referencia

r=filter((1-af),[1 -af],w); %os dados do vetor w são filtrados por um filtro de 1 ordem
% F(z)= (1-af)/(z-af), descrito pelo numerador (1-af) e denominador [1 -af].

% inicializa

inc_u=0;


% começo do laço de controle

for k=1:pontos
   
% calculo da saída do processo

    yp(k)=Ap*y_ant'+Bp*u_ant'; % y processo

% perturbação determinística para a saída
if k>150 % a pert entra entre os pontos 150 e 220
 ypr(k)=yp(k)+0.1; %perturbação constante tipo degrau
else
 ypr(k)=yp(k);  %ypr é a variavel de processo com a perturbacao
end


% atualiza a saída para simulação do proximo passo

if na==1
    y_ant=yp(k);
else
    aux_y=y_ant(1:na-1); % aux_y recebe do primeiro ao penultimo elemento de y_ant
    y_ant=[yp(k) aux_y]; % y_ant recebe o yp(atual) e da um shift de segundo elemento até o ultimo com os elementos passados de y_ant( descarta o ultimo elemento )
end




% Calcula a predicao para o DMC

   f=zeros(1,p); % inicializa o vetor f com 'p' zeros. (20 no caso)

      for kk=1:p
% Monta um vetor com as parcelas (gkk+i - gkk)
% o calculo se faz em dois passos porque quando o indice kk+i passa de N,
% nao teriamos valores no vetor g para esse valor, mas como se assume que a planta é estavel
% e N é grande, todos os g(i) = g(N) para i>N

      for i=1:N-kk % N periodos de amostragem que a planta se estabiliza, p é o horizonte
             vect_g(i)=g(kk+i)-g(i); 
         end
         for i=N-kk+1:N
             vect_g(i)=g(N)-g(i);
         end
         f(kk)=ypr(k)+vect_g*u_livre';
      end

% considera ref futura cte

ref=r(k)*ones(1,p)';


% calcula o controle

inc_u=qn*(ref-f'); % qn = k1, deltaU = k(W - f)
if k==1
    u(k)=inc_u;
else
    u(k)=u(k-1)+ inc_u;
end

% atualizacao do vetor de controle

aux_u=u_ant(1:nb-1);
u_ant=[u(k) aux_u];

%atualizacao de  u_livre
% u_livre= [du(k-1) du(k-2) ..... du(k-N)]

aux_2=u_livre(1:N-1);
u_livre=[inc_u aux_2];


end


nm=pontos;

% parameters for figures
tamletra=12;
tamnum=12; 

h=subplot(2,1,1);
plot(ypr(1:nm),'-');
hold
%plot(t(1:nm),y1(1:nm),'-');
plot(w(1:nm),'--');


ylabel('output', 'FontSize', tamletra);
xlabel('time', 'FontSize', tamletra);
%title('(b)', 'FontSize', tamletra);
%axis([0 pontos -0.1 0.7])


h=subplot(2,1,2);
plot(u(1:nm),'-');
hold
%plot(t(1:nm),u1(1:nm),'-');
ylabel('control action', 'FontSize', tamletra);
xlabel('time', 'FontSize', tamletra);
%axis([0 pontos -2 1])
