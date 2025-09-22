
% Script MATLAB para filtro LCL passivamente amortecido
% LCL com elemento resisitivo alocado em paralelo com o capacitor de
% filtragem

L1 = 0.001;
Lc = 0.0005;
L2 = 0.001;
Rf = 0.1;
Rd = 10;
Cf = 1e-05;

omega_d = sqrt((L1 + Lc)/(L2 * Rd * Cf));
z_d = 1/(2 * omega_d * Rd * Cf);

num = [omega_d^2];
den = [1, 2*z_d*omega_d, omega_d^2];
GvL = tf(num, den);
bode(GvL);
grid on;


% Parâmetros fornecidos
Le = 0.52e-3;   % Indutância equivalente (H)
Cf = 3.9e-6;    % Capacitância do filtro (F)
Vuf = 220;      % Tensão aplicada (V)

% Cálculo da resistência de amortecimento Rd
Rd = 0.5 * sqrt(Le / Cf);

% Corrente através do resistor
Ird = Vuf / Rd;

% Potência dissipada no resistor
Prd = Vuf^2 / Rd;

% Exibir resultados
fprintf('Rd = %.2f Ohms\\n', Rd);
fprintf('Ird = %.2f A\\n', Ird);
fprintf('Prd = %.2f kW\\n', Prd / 1000);

% Função de transferência para análise de Bode
s = tf('s');
GvL = (Rd * Cf * s^2 + (1 + Rd/Rd) * Cf * s + 1)^(-1);  % Simplificação da equação (5.50)

% Plot da resposta em frequência
figure;
bode(GvL);
grid on;
title('Resposta em Frequência - Filtro LCL com Amortecimento Passivo');


%% 

% Script MATLAB para filtro LCL passivamente amortecido
% LCL com elemento resisitivo alocado em paralelo com o capacitor de
% filtragem
