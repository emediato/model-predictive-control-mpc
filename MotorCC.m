clear all
close all
clc

% Graficos 
t = (0:0.001:2.5);   % resultados no intervalo [0, 2.5] com resolucao 0.001

x = (t>=0);              % degrau unitario

y0 = -0.0617*exp(-115.66*t) + 1.0617*exp(-4.84*t);
yesn = 0.9318*(1-exp(-4.84*t)) - 0.0390*(1-exp(-115.66*t));
y = y0 + yesn;

figure
plot(t,x,'g',t,y0,'r',t,yesn,'b',t,y,'k');
legend('x(t)','y0(t)','yesn(t)','y(t)=y0(t)+yesn(t)');
grid on;

figure
plot(t,x,'g',t,y,'k');
legend('x(t)','y(t)=y0(t)+yesn(t)');
grid on;


% Simulacao no Espaco de Estados
v0 = [1 2];              % condicao inicial do Sistema Original

A = [-120.5 1; -560 0];
B = [0 ; 500];
C = [1 0];
D = [0];

K = ctrb(A',C')';

q0 = inv(K)*v0';     % condicao inicial no Espaco de Estados
