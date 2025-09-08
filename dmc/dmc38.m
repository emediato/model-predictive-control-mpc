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

begin_g = 0.1;
last_g = 1;

ts = 1;

rho = 1;
delta = 0.5;

qtt = 4;
g0 = linspace(begin_g, last_g, qtt);

g0=[0.1000,    0.3000,    0.7000, 0.9];
g0 = [g0,1,1,1,1,1];

% u(k) = alfa1 * u (k-1) + alfa2 + ... + beta * y(k) + gamma * r(k);
coef_g = 1/delta * g0; 
g=coef_g;

% horizonte controle
Nu = 2 ;
% horizonte predicao
N = 3;
 
Nss = 5;

G1 = coef_g(1:N);
G2 = coef_g(1:Nu);

G = [G1;0,G2];
G = G';

size_eye = length(G) ;

% GAIN
K = G'*G ;

size_eye = length(K) ;
%%
% K dimension is [Nu, N]

K = K + delta*eye(size_eye) ;
K = inv(K) * G' ;

K1=K(1,1:end);

k1 = K1(1); k2 = K1(2) ; k3=K1(3);

% primeira linha de K
disp(G'*G + 0.5);


f = zeros(3,1);

syms k
syms r(k) y(k) 
syms u_k_minus_1 u_k_minus_2 u_k_minus_3 u_k_minus_4 u_k_minus_5
syms u_k delta_u_k % current control


%% lei de controle
% delta_u_k = k1 * (r(k) - f(k+1)) + k2 * (r(k) - f(k+2)) + k3 * (r(k) - f(k+3));

% resposta livre para os 3 elementos do horizonte de predicao
% f(k+1) = (g(2)-g(1))*delta_u(k-1) + (g(3)-g(2))*delta_u(k-2) + ...
%             + (g(4)-g(3))*delta_u(k-3) + (g(5)-g(4))*delta_u(k-4) + ...
%                                             (g(6)-g(5))*delta_u(k-5) + y(k)
% 
% f(k+2) = (g(3)-g(1))*delta_u(k-1) + (g(4)-g(2))*delta_u(k-2) + ...
%             + (g(5)-g(3))*delta_u(k-3) + (g(7)-g(4))*delta_u(k-4) + ...
%                                             (g(8)-g(5))*delta_u(k-5) + y(k)
% 
% f(k+3) = (g(4)-g(1))*delta_u(k-1) + (g(5)-g(2))*delta_u(k-2) + ...
%             + (g(6)-g(3))*delta_u(k-3) + (g(7)-g(4))*delta_u(k-4) + ...
%                                             (g(8)-g(5))*delta_u(k-5) + y(k)
f_k_1 = (g(2)-g(1))*u_k_minus_1 + (g(3)-g(2))*u_k_minus_2+ ...
            + (g(4)-g(3))*u_k_minus_3 + (g(5)-g(4))*u_k_minus_4 + ...
                                            (g(6)-g(5))*u_k_minus_5 + y(k)

f_k_2 = (g(3)-g(1))*u_k_minus_1 + (g(4)-g(2))*u_k_minus_2 + ...
            + (g(5)-g(3))*u_k_minus_3 + (g(6)-g(4))*u_k_minus_4 + ...
                                            (g(7)-g(5))*u_k_minus_5 + y(k)

f_k_3 = (g(4)-g(1))*u_k_minus_1 + (g(5)-g(2))*u_k_minus_2+ ...
            + (g(6)-g(3))*u_k_minus_3 + (g(7)-g(4))*u_k_minus_4 + ...
                                            (g(8)-g(5))*u_k_minus_5 + y(k)

% lei de controle
%delta_u(k) = k1 * (r(k) - f(k+1)) + k2 * (r(k) - f(k+2)) + k3 * (r(k) - f(k+3));


delta_u_k = k1 * (r(k) - y(k)) + k2 * (r(k) - y(k)) + k3 * (r(k) - y(k))... 
   - k1 * ( f_k_1) - k2 * (f_k_2) - k3 * ( f_k_3)


% collect(simplify(delta_u_k), u_k_minus_1)



% Coefficients calculated with high precision
coeff_u_k_minus_1   = 440/459;      % ≈ 0.9586056644880175
coeff_u_k_minus_2   = 422/459;      % ≈ 0.9193899773420471
coeff_u_k_minus_3   = 64/153;       % ≈ 0.4183006535947713
coeff_u_k_minus_4   = 202/1377;     % ≈ 0.1466957153231663
y_k_coeff     = -2020/1377;    % ≈ 1.4669571532316630  % Note: 2020 instead of 1010
r_k_coeff     = 1010/1377;    % ≈ 0.7334785766158315

coeff_r_minus_k= r_k_coeff+y_k_coeff;
%G_U_E =collect(simplify(u_k, r(k)) 



%%
% delta_u(k) = u(k) - u(k-1)
% u(k) = delta_u(k) + u(k-1)

% u_k = delta_u_k + u_k_minus_1 ;


u_k = u_k_minus_1 + coeff_r_minus_k * (r(k) - y(k)) ...
        - coeff_u_k_minus_1 *( u_k_minus_1 -  u_k_minus_2) ...
             - coeff_u_k_minus_2 *( u_k_minus_2 -  u_k_minus_3) ...
                - coeff_u_k_minus_3 *( u_k_minus_3 -  u_k_minus_4) ...
                    - coeff_u_k_minus_4 *( u_k_minus_4 -  u_k_minus_5)
