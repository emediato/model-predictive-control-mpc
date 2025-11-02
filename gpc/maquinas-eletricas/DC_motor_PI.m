% Parameters from the first image
V_in = 48;          % Bus voltage [V]
R_a = 0.016;        % Armature resistance [Ω]
v_t = 18;            % armature voltage
i_a = 100;           % armature current [A]
i_f = 1.6 ;           %  field current [A]
V_trip = 3.3;         % Triangular peak value [V]
L_a = 0.06; % indutancia de armadura


% Armature mutual inductance [H]
L_af = (v_t - R_a*i_a )/ 70 * v_t * 2*pi/60 * i_f;


% Controller gains
k_i = 1;            % Current sensor gain
kia = 1;            % Current controller gain (?)
kpwm = 1/V_trip;    % PWM gain

% Transfer function G_id(s) = i_a(s)/d(s) = V_in / (R_a + s*L_a)
% Create the transfer function
numerator = V_in;
denominator = [L_a, R_a];
G_id = tf(numerator, denominator);

% Open-loop transfer function FTLA
FTLA_net = k_i * kpwm * G_id;

% Create frequency vector 
fs = 1.66e4;           % Frequência do conversor [Hz]
fci = fs / 10;         % Frequência de cruzamento [Hz]
 w_cross = 2 * pi * fci;     % Frequência angular de cruzamento [rad/s]

omega =  w_cross;

% Calculate magnitude and phase
[mag, phase] = bode(FTLA_net, omega);
fprintf('Na frequência ωc = %.2f rad/s:\n', omega);
fprintf('  |FTLA(ω)| = %.4f\n', mag);
fprintf('  ∠FTLA(ω) = %.2f°\n', phase);

phase_NC = (phase)* 180/pi;
fprintf('  ∠FTLA(ω)| = %.4f \n', phase_NC);

fprintf('Frequência de cruzamento:  fci = %.2f Hz\n',  fci);
fprintf('Frequência angular de cruzamento:  w_cross = %.2f rad/s\n',  w_cross);

% Margem de fase desejada
M_p = 60 * pi/180;      % 60 graus em radianos
fprintf('Margem de fase desejada: M_p = %.2f graus, %.2f rad/s\n', M_p*180/pi, M_p);

% Função de transferência da planta de corrente
G_id = tf(V_in, [L_a, R_a]);
fprintf('\nFunção de transferência G_id(s):\n');
disp(G_id);

% FTLA inicial (malha aberta sem compensador)
FTLA_int = k_i * kpwm * G_id;
fprintf('FTLA inicial (sem compensador):\n');
disp(FTLA_int);

% Avaliar FTLA_int na frequência  w_cross
[mag_wd, phase_wd] = bode(FTLA_int,  w_cross);

fprintf('|FTLA_int( w_cross)| = %.4f\n', mag_wd);
fprintf('arg(FTLA_int( w_cross)) = %.2f graus\n', phase_wd);

sys=feedback(FTLA_int,1);
[mag_wd, phase_wd] = bode(sys,  w_cross);


fprintf('\nNa frequência  w_cross = %.2f rad/s:\n',  w_cross);
fprintf('|FTLA_int( w_cross)| = %.4f\n', mag_wd);
fprintf('arg(FTLA_int( w_cross)) = %.2f graus\n', phase_wd);
angPI = (M_p - pi - phase_wd*pi/180)*180/pi ;


% Cálculo do zero do compensador 
w_z =  w_cross/tan(M_p - pi/2 - phase_wd*pi/180); % Fase necessária do compensador
tau_z = 1/w_z;
fprintf('Frequência do zero do compensador: w_z = %.2f rad/s\n', w_z);

k_ci =  w_cross/((sqrt(w_z^2 +  w_cross^2)*(mag)))

% Controlador PI: C_i(s) = k_ci * (s + w_z) / s
C_i = k_ci * tf([1, w_z], [1, 0]);

fprintf('\nControlador PI C_i(s):\n');
disp(C_i);

% FTLA com compensador
FTLA_comp = FTLA_int * C_i;

% Malha fechada de corrente
G_cl_current = feedback(FTLA_comp, 1);
fprintf('Função de transferência de malha fechada de corrente:\n');
disp(G_cl_current);


% 
% % Create Bode plot
% figure;
% 
% % Magnitude plot
% subplot(2,1,1);
% semilogx(omega, mag_dB, 'b', 'LineWidth', 2);
% grid on;
% xlabel('Frequency (rad/s)');
% ylabel('Magnitude (dB)');
% title('Bode Diagram - Magnitude');
% legend('FTLA_{net}');
% 
% % Phase plot
% subplot(2,1,2);
% semilogx(omega, phase, 'r', 'LineWidth', 2);
% grid on;
% xlabel('Frequency (rad/s)');
% ylabel('Phase (degrees)');
% title('Bode Diagram - Phase');
% legend('FTLA_{net}');
% 
% % Alternative: Use MATLAB's built-in bode plot
% figure;
% bode(FTLA_net);
% grid on;
% title('Bode Diagram using MATLAB bode function');
% 
% % Display transfer function information
% disp('Transfer Function G_id(s):');
% disp(G_id);
% disp('Open-loop Transfer Function FTLA_net(s):');
% disp(FTLA_net);
% disp(['DC Gain: ', num2str(dcgain(FTLA_net))]);
% disp(['Crossover Frequency: ', num2str(bandwidth(FTLA_net)), ' rad/s']);
% 
% 

%%  MALHA DE VELOCIDADE
fprintf('\n=== MALHA DE VELOCIDADE ===\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = 0.018;              % Momento de inércia [kg·m²]
fcv = 1.66e3/100;
w_cv = 2 * pi * fcv;     % Frequência angular de cruzamento [rad/s]

k_i = 1;                % sensor de corrente
kia = 1;                % controlador de corrente
kpwm = 1;        % Ganho PWM
% Margem de fase desejada
M_p = 60 * pi/180;      % 60 graus em radianos

% Função de transferência da planta de velocidade
% G_n(w) = (30/pi) * (L_af * I_f) / (j * w * J)
% Convertendo para função de transferência no domínio s:
% G_n(s) = K / s, onde K = (30/pi) * (L_af * I_f) / J
K = ((30/pi) * L_a * i_f)/J;
num = K;
den = [1 0];
G_n =  tf( [num, den]);
fprintf('\nFunção de transferência G_n(s):\n');
disp(G_n);

K_vel = (30/pi) * (L_a * i_f) / J;
fprintf('Ganho da planta de velocidade: K_vel = %.4f\n', K_vel);

G_vel = tf(K_vel, [1, 0]);
fprintf('Função de transferência G_vel(s):\n');
disp(G_vel);

% Considerando k_s = 1 (ganho do sensor de velocidade)
k_s = 1;
FILA_vel = (1/k_s) * k_s * G_vel;  % FTLA para malha de velocidade
fprintf('FTLA para malha de velocidade:\n');
disp(FILA_vel);

% FTLA inicial (malha aberta sem compensador)
FTLA_NC_vel = k_i * kpwm * G_n;

% Avaliar FTLA_int na frequência  w_cross
[mag_vel, phase_vel] = bode(FILA_vel, w_cv);

angPIvel = (M_p - pi - phase_vel*pi/180)*180/pi 

% Cálculo do zero do compensador 
w_v = w_cv/(tan(M_p - pi/2 - phase_vel * pi/180)); % Fase necessária do compensador

tau_cv = 1/w_v;
fprintf('Frequência do zero do compensador: w_v = %.2f rad/s\n', w_v);

% Cálculo do ganho k_d
k_cv = w_cv / (sqrt(w_v^2 + w_cv^2)*mag_vel);
fprintf('Ganho do compensador: k_d = %.4f\n', k_cv);


