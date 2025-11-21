
% Dados do sistema
fs = 20000;           % Frequência de chaveamento (Hz)
fciq = fs/10;         % Frequência de corte (Hz)
wciq = 2*pi*fciq;     % Frequência angular de corte (rad/s)
Mphi = pi/3;          % Margem de fase (rad)

% Função de transferência não compensada
FTLA_nc = tf(8333.33, [1 100]);  % FTLA não compensada

% Cálculo do zero do compensador
% Aproximação: ângulo de FTLA_nc em wciq
[mag, phase] = bode(FTLA_nc, wciq);
phase = deg2rad(phase);  % Converter para radianos
wziq = wciq / tan(Mphi - pi/2 - phase);

% Ganho do compensador
Kciq = (wciq / sqrt(wciq^2 + wziq^2)) * (1/mag);

% Compensador
Ciq = Kciq * tf([1 wziq], [1 0]);

% FTLA compensada
FTLA_c = series(Ciq, FTLA_nc);

% Exibir resultados
disp('Zero do compensador (rad/s):'); disp(wziq);
disp('Ganho do compensador:'); disp(Kciq);
disp('Compensador Ciq(s):'); Ciq
disp('FTLA compensada:'); FTLA_c

% Plot de Bode
figure;
bode(FTLA_nc, FTLA_c);
legend('Não compensada', 'Compensada');
grid on;


%% new
% Dados
fs = 20000;                % Frequência de chaveamento (Hz)
fciq = fs/10;              % Frequência de corte (Hz)
wciq = 2*pi*fciq;          % Velocidade angular de corte (rad/s)
Mphi = pi/3;               % Margem de fase (rad)

% Função de transferência não compensada
FTLA_nc = tf(8333.33, [1 100]);

% Ângulo da FTLA_nc em wciq
[mag, phase] = bode(FTLA_nc, wciq);
phase = deg2rad(phase);    % Converter para radianos

% Zero do compensador
wziq = wciq / tan(Mphi - pi/2 - phase);

% Ganho do compensador
Kciq = (wciq / sqrt(wciq^2 + wziq^2)) * (1/mag);

% Compensador
Ciq = Kciq * tf([1 wziq], [1 0]);

% FTLA compensada
FTLA_c = series(Ciq, FTLA_nc);

% Exibir resultados
disp('Zero do compensador (rad/s):'); disp(wziq);
disp('Ganho do compensador:'); disp(Kciq);
disp('Compensador Ciq(s):'); Ciq
disp('FTLA compensada:'); FTLA_c

% Plot de Bode
figure;
bode(FTLA_nc, FTLA_c);
legend('Não compensada', 'Compensada');
grid on;
