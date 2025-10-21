%% Controle Preditivo com Polinômio C - Exemplo Prático
clear; clc; close all;

%% Definição do sistema
A = [1, -1.2, 0.35];    % A(z^-1) = 1 - 1.2z^-1 + 0.35z^-2
B = [0.4, 0.1];         % B(z^-1) = 0.4 + 0.1z^-1  
C = [1, 0.5, -0.2];     % C(z^-1) = 1 + 0.5z^-1 - 0.2z^-2 (MÔNICO)
d = 1;                  % Atraso

%% Parâmetros de predição
N1 = 1;     % Horizonte inferior
N2 = 3;     % Horizonte superior
Nu = 2;     % Horizonte de controle

%% Simulação
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
    xi(k) = -C(2)*xi(k-1) - C(3)*xi(k-2) + e(k);
end

% Saída do sistema
y = zeros(Nsim, 1);
y_clean = zeros(Nsim, 1);  % Sem ruído

for k = 3:Nsim
    % Sem ruído
    y_clean(k) = -A(2)*y_clean(k-1) - A(3)*y_clean(k-2) + ...
                 B(1)*u(k-d) + B(2)*u(k-d-1);
    
    % Com ruído
    y(k) = -A(2)*y(k-1) - A(3)*y(k-2) + ...
            B(1)*u(k-d) + B(2)*u(k-d-1) + xi(k);
end

%% Cálculo dos polinômios E_j e F_j considerando C
% Para j = 1
na = length(A) - 1;
nb = length(B) - 1;
nc = length(C) - 1;

% Polinômio A tilde
A_tilde = conv([1 -1], A);  % ΔA = (1-z^-1)A(z^-1)

% Resolver equação diofantina: C = E_j * A_tilde + z^-j * F_j
% Para j = 1
E1 = [1];  % E1 tem grau 0
% C = E1*A_tilde + z^-1*F1
% F1 = z*(C - E1*A_tilde)

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

%%
C_minus_E1A = C - conv(E1, A_tilde);
% Remover o primeiro elemento (coeficiente de z^0) para obter F1
F1 = C_minus_E1A(2:end);

fprintf('=== Verificação Equação Diofantina j=1 ===\n');
fprintf('E1 = [%.4f]\n', E1);
fprintf('F1 = [%.4f', F1);
if length(F1) > 1
    for i = 2:length(F1)
        fprintf(' %.4f', F1(i));
    end
end
fprintf(']\n');

% Verificação
lado_esquerdo = C;
lado_direito = conv(E1, A_tilde) + [0, F1];
fprintf('C:                [%.4f %.4f %.4f]\n', C);
fprintf('E1*A_tilde + F1: [%.4f %.4f %.4f]\n', lado_direito);
fprintf('Erro: %.6f\n', norm(lado_esquerdo - lado_direito));

%% Predição um passo à frente usando C
y_pred = zeros(Nsim, 1);

for k = 4:Nsim
    % Predição ŷ(k|k-1) = [E1*B/C] Δu(k-1) + [F1/C] y(k-1)
    
    % Numerador E1*B
    num = conv(E1, B);
    
    % Usar filtro recursivo para (E1*B)/C e F1/C
    % Simplificação: aproximação por resposta ao degrau
    delta_u = [u(k-d) - u(k-d-1)];
    
    % Parte forçada (simplificada)
    forced_part = filter(num, C, delta_u);
    
    % Parte livre
    free_part = filter(F1, C, y(k-1:-1:max(1,k-3)));
    
    y_pred(k) = forced_part(1) + free_part(1);
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

%% Função para equação diofantina generalizada
fprintf('\n=== Cálculo dos Polinômios E_j e F_j ===\n');
for j = 1:3
    [Ej, Fj] = diophantine_eq(A_tilde, C, j);
    fprintf('j = %d:\n', j);
    fprintf('  E%d = [', j); fprintf('%.4f ', Ej); fprintf(']\n');
    fprintf('  F%d = [', j); fprintf('%.4f ', Fj); fprintf(']\n');
    
    % Verificação
    lado_esq = C;
    lado_dir = conv(Ej, A_tilde);
    lado_dir(end-j+1:end) = lado_dir(end-j+1:end) + Fj;
    erro = norm(lado_esq - lado_dir);
    fprintf('  Erro verificação: %.6f\n', erro);
end

function [Ej, Fj] = diophantine_eq(A_tilde, C, j)
    % Resolve C = Ej * A_tilde + z^-j * Fj
    % Ej tem grau j-1, Fj tem grau grau(A_tilde)-1
    
    nA = length(A_tilde) - 1;  % grau de A_tilde
    nC = length(C) - 1;        % grau de C
    
    % Grau de Ej = j-1
    % Grau de Fj = nA-1
    
    % Montar sistema de equações
    nEq = nC + 1;
    M = zeros(nEq, j + nA);
    
    % Preencher matriz para Ej * A_tilde
    for i = 1:j
        for m = 1:length(A_tilde)
            if (i + m - 1) <= nEq
                M(i + m - 1, i) = M(i + m - 1, i) + A_tilde(m);
            end
        end
    end
    
    % Preencher matriz para z^-j * Fj
    for i = 1:nA
        for m = 0:0
            if (j + i + m - 1) <= nEq
                M(j + i + m, j + i) = M(j + i + m, j + i) + 1;
            end
        end
    end
    
    % Resolver sistema
    b = C(:);
    x = M \ b;
    
    Ej = x(1:j)';
    Fj = x(j+1:end)';
end