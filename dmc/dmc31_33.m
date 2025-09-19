% Codigo que implementa DMC

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




 %%
len = Nu;

delta_u = repmat([0]',len-1,1);
h_impulse = repmat([0]',len-1,1);

for k = 1:len-1
    delta_u(k) = y(k+1) - g(k);

    % resposta ao impulso
    % h(k) = g(k) - g(k-1);
    %h(k) = g(k) - g(k-1)
    h_impulse(k+1) = g(k+1) - g(k); % delta_g


    % resposta do sistema - pagina 47
    y_new = (h_impulse(k) * delta_u(k))/G;
end
