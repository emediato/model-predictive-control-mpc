%% Controle Preditivo com Polinômio C - Exemplo Prático
clear; clc; close all;

%% Definição do sistema
A = [1, -.9];    
B = [0.3];         
C = [1, 0.4];     % C(z^-1) = 1 + 0.4z^-1  (MÔNICO)
d = 1;                  % Atraso

%% Parâmetros de predição
N1 = 1;     % Horizonte inferior
N2 = 3;     % Horizonte superior
Nu = 2;     % Horizonte de controle

%% Cálculo dos polinômios E_j e F_j considerando C
fprintf('=== Cálculo dos Polinômios E_j e F_j ===\n');

%%

% Polinômio A tilde
A_tilde = conv([1 -1], A);  % ΔA = (1-z^-1)A(z^-1)
fprintf('A_tilde = [');
fprintf('%.4f ', A_tilde);
fprintf(']\n');

% Para j = 1
E1 = [1];  % E1 tem grau 0

% C = E1*A_tilde + z^-1*F1
% Precisamos igualar os tamanhos dos vetores
E1A = conv(E1, A_tilde);

% Ajustar tamanhos para subtração
max_len = max(length(C), length(E1A));
C_padded = [C, zeros(1, max_len - length(C))];
E1A_padded = [E1A, zeros(1, max_len - length(E1A))];

C_minus_E1A = C_padded - E1A_padded;

% F1 são os coeficientes de z^-1 em diante
F1 = C_minus_E1A(2:end);

fprintf('\n--- Para j = 1 ---\n');
fprintf('E1 = [%.4f]\n', E1);
fprintf('F1 = [');
fprintf('%.4f ', F1);
fprintf(']\n');

% Verificação
lado_esquerdo = C;
lado_direito = conv(E1, A_tilde);
lado_direito = [lado_direito, zeros(1, max_len - length(lado_direito))];
lado_direito(2:end) = lado_direito(2:end) + F1;

% Remover zeros à direita
lado_direito = lado_direito(1:length(C));

fprintf('C:                [');
fprintf('%.4f ', C);
fprintf(']\n');
fprintf('E1*A_tilde + F1: [');
fprintf('%.4f ', lado_direito);
fprintf(']\n');
fprintf('Erro: %.6f\n', norm(C - lado_direito));




% Para j = 1
E1 = [1];
F1 = C(2) - A(2);  % F1 = -0.4 - (-0.9) = 0.5
fprintf('j = 1:\n');
fprintf('  E1 = [%.4f]\n', E1);
fprintf('  F1 = [%.4f]\n', F1);

% Verificação j=1
verif1 = conv(E1, A);
verif1(2) = verif1(2) + F1;
fprintf('  Verificação: C = [%.4f %.4f], E1*A + F1 = [%.4f %.4f]\n', ...
        C(1), C(2), verif1(1), verif1(2))


%% Para j = 2
fprintf('\n--- Para j = 2 ---\n');

% E2 tem grau 1: E2 = e20 + e21*z^-1
% C = E2*A_tilde + z^-2*F2

% Montar sistema de equações:
% Coeficientes de z^0: 1 = e20 * 1
% Coeficientes de z^-1: c1 = e20*a1 + e21*1
% Coeficientes de z^-2: c2 = e20*a2 + e21*a1 + f20
% Coeficientes de z^-3: c3 = e20*a3 + e21*a2 + f21

e20 = 1;  % de z^0

% De z^-1: c1 = e20*a1 + e21
e21 = C(2) - e20 * A_tilde(2);

% De z^-2: c2 = e20*a2 + e21*a1 + f20
f20 = C(2) - e20 * A_tilde(3) - e21 * A_tilde(2);

% De z^-3: c3 = e20*a3 + e21*a2 + f21
f21 = 0 - e20 * 0 - e21 * A_tilde(3);  % c3 = 0 para nosso C

E2 = [e20, e21];
F2 = [f20, f21];

fprintf('E2 = [');
fprintf('%.4f ', E2);
fprintf(']\n');
fprintf('F2 = [');
fprintf('%.4f ', F2);
fprintf(']\n');

% Verificação
lado_dir2 = conv(E2, A_tilde);
lado_dir2 = [lado_dir2, zeros(1, length(C) - length(lado_dir2))];
lado_dir2(3:end) = lado_dir2(3:end) + F2;
fprintf('E2*A_tilde + F2: [');
fprintf('%.4f ', lado_dir2);
fprintf(']\n');
fprintf('Erro: %.6f\n', norm(C - lado_dir2(1:length(C))));

%% Simulação do sistema
fprintf('\n=== Simulação do Sistema ===\n');

Nsim = 100;
t = 1:Nsim;

% Entrada (degraus)
u = zeros(Nsim, 1);
u(20:40) = 1;
u(60:80) = -0.5;

% Ruído colorido (filtrado por C)
e = 0.1 * randn(Nsim, 1);  % Ruído branco
xi = zeros(Nsim, 1);       % Ruído colorido

for k = 3:Nsim
    if k > 2
        xi(k) = e(k) + C(2)*xi(k-1) ; %+ C(3)*xi(k-2);
    end
end

% Saída do sistema
y = zeros(Nsim, 1);
y_clean = zeros(Nsim, 1);  % Sem ruído

for k = 3:Nsim
    % Sem ruído
    y_clean(k) = -A(2)*y_clean(k-1) + ... % - A(3)*y_clean(k-2) + ...
                 B(1)*u(max(1,k-d)) ; %+ B(2)*u(max(1,k-d-1));
    
    % Com ruído
    y(k) = -A(2)*y(k-1) + ... - A(3)*y(k-2) + ...
            B(1)*u(max(1,k-d))  + xi(k); % + B(2)*u(max(1,k-d-1))
end

%% Predição um passo à frente
y_pred = zeros(Nsim, 1);

for k = 4:Nsim
    % Forma incremental do modelo:
    % A_tilde * y(k) = B * Δu(k-1) + C * e(k)
    
    % Predição: ŷ(k|k-1) = [B/A_tilde] Δu(k-1) + [1 - C/A_tilde] y(k-1)
    
    % Parte forçada
    delta_u = u(k-d) - u(k-d-1);
    forced_part = filter(B, A_tilde, delta_u);
    
    % Parte livre (aproximação)
    free_part = y(k-1) - filter(C, A_tilde, y(k-1));
    
    y_pred(k) = forced_part + free_part;
end

%% Gráficos
figure('Position', [100 100 1200 800]);

subplot(3,1,1);
plot(t, u, 'b-', 'LineWidth', 2);
title('Sinal de Controle u(k)');
xlabel('Amostras');
ylabel('u(k)');
grid on;

subplot(3,1,2);
plot(t, y_clean, 'g--', 'LineWidth', 1.5); hold on;
plot(t, y, 'r-', 'LineWidth', 1.5);
plot(t, y_pred, 'b:', 'LineWidth', 1.5);
legend('Sem Ruído', 'Com Ruído', 'Predição 1-passo', 'Location', 'best');
title('Saída do Sistema e Predição');
xlabel('Amostras');
ylabel('y(k)');
grid on;

subplot(3,1,3);
plot(t, xi, 'm-', 'LineWidth', 1.5);
title('Ruído Colorido ξ(k) = C(z^{-1})e(k)');
xlabel('Amostras');
ylabel('ξ(k)');
grid on;

% Testar a função
% for j = 1:3
%     [Ej, Fj] = diophantine_eq(A_tilde, C, j);
%     fprintf('j = %d:\n', j);
%     fprintf('  E%d = [', j); fprintf('%.4f ', Ej); fprintf(']\n');
%     fprintf('  F%d = [', j); fprintf('%.4f ', Fj); fprintf(']\n');
% 
%     % Verificação
%     lado_dir = conv(Ej, A_tilde);
%     lado_dir_padded = [lado_dir, zeros(1, length(C) - length(lado_dir))];
%     lado_dir_padded(j+1:j+length(Fj)) = lado_dir_padded(j+1:j+length(Fj)) + Fj;
%     erro = norm(C - lado_dir_padded(1:length(C)));
%     fprintf('  Erro: %.6f\n', erro);
% end

%% Função melhorada para equação diofantina
fprintf('\n=== Cálculo com Função Diofantina ===\n');

function [Ej, Fj] = diophantine_eq(A_tilde, C, j)
    % Resolve C = Ej * A_tilde + z^-j * Fj
    
    nA = length(A_tilde) - 1;  % grau de A_tilde
    nC = length(C) - 1;        % grau de C
    
    % Número total de coeficientes desconhecidos
    n_unknowns = j + nA;  % j coeficientes em Ej + nA coeficientes em Fj
    
    % Montar sistema de equações
    M = zeros(nC + 1, n_unknowns);
    b = C(:);
    
    % Preencher matriz para Ej * A_tilde
    for i = 1:j  % coeficientes de Ej
        for m = 1:length(A_tilde)
            row = i + m - 1;
            if row <= nC + 1
                M(row, i) = M(row, i) + A_tilde(m);
            end
        end
    end
    
    % Preencher matriz para z^-j * Fj
    for i = 1:nA  % coeficientes de Fj
        for m = 1:1  % Fj é multiplicado por z^-j
            row = j + i + m - 1;
            if row <= nC + 1
                M(row, j + i) = M(row, j + i) + 1;
            end
        end
    end
    
    % Resolver sistema
    x = M \ b;
    
    Ej = x(1:j)';
    Fj = x(j+1:end)';
end
% 
