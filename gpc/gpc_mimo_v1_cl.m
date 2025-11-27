%% ========================================================================
%  MODELO CORRETO: GPC MIMO PARA CONVERSOR TRIFÁSICO
%  Entradas: [Vd; Vq] | Saídas: [Id; Iq]
%% ========================================================================

clear all; close all; clc;

%% PARÂMETROS DO HARDWARE
%==========================================================================
Ts = 1/1000;       % Período de amostragem [s] - 1 kHz
Vdc = 60;          % Tensão do barramento CC [V]
R = 22.5;          % Resistência de fase [Ω]
L = 0.10;          % Indutância de fase [H]
x_ref = [0.5; 0.5];% Referências [Id_ref; Iq_ref] [A]

% Frequência da rede
f_rede = 60;       % Hz
omega_e = 2*pi*f_rede; % rad/s

fprintf('=== PARÂMETROS DO SISTEMA ===\n');
fprintf('R = %.2f Ω, L = %.4f H\n', R, L);
fprintf('τ = L/R = %.4f s (%.2f ms)\n', L/R, L/R*1000);
fprintf('ωe = %.2f rad/s (%.1f Hz)\n', omega_e, f_rede);
fprintf('Vdc = %.1f V\n', Vdc);

%% MODELO EM ESPAÇO DE ESTADOS
%==========================================================================
I_2 = eye(2);

%--- OPÇÃO 1: MODELO DESACOPLADO (Simplificado) ---
fprintf('\n=== OPÇÃO 1: MODELO DESACOPLADO ===\n');
fprintf('Considera ωe constante, acoplamento tratado como perturbação\n\n');

% Matrizes contínuas
F_desacoplado = -(R/L) * I_2;
G_desacoplado = (1/L) * I_2;

fprintf('Matriz F (desacoplado):\n');
disp(F_desacoplado);
fprintf('Matriz G (desacoplado):\n');
disp(G_desacoplado);

% Discretização
A_desacoplado = expm(F_desacoplado * Ts);
B_desacoplado = -F_desacoplado^(-1) * (I_2 - A_desacoplado) * G_desacoplado;

fprintf('Matriz A discreta:\n');
disp(A_desacoplado);
fprintf('Matriz B discreta:\n');
disp(B_desacoplado);

%--- OPÇÃO 2: MODELO ACOPLADO (Completo) ---
fprintf('\n=== OPÇÃO 2: MODELO ACOPLADO ===\n');
fprintf('Inclui acoplamento cruzado dq explicitamente\n\n');

% Matrizes contínuas com acoplamento
F_acoplado = [-R/L,     omega_e;
              -omega_e, -R/L];
          
G_acoplado = [1/L,  0;
              0,    1/L];

fprintf('Matriz F (acoplado):\n');
disp(F_acoplado);
fprintf('Matriz G (acoplado):\n');
disp(G_acoplado);

% Discretização
A_acoplado = expm(F_acoplado * Ts);
B_acoplado = F_acoplado^(-1) * (A_acoplado - I_2) * G_acoplado;

fprintf('Matriz A discreta (acoplado):\n');
disp(A_acoplado);
fprintf('Matriz B discreta (acoplado):\n');
disp(B_acoplado);

%% FUNÇÕES DE TRANSFERÊNCIA G(i,j)
%==========================================================================
fprintf('\n=== FUNÇÕES DE TRANSFERÊNCIA ===\n');

% Definição de variáveis
m = 2; % número de ENTRADAS: [Vd, Vq]
n = 2; % número de SAÍDAS: [Id, Iq]
mq = 2; % número de PERTURBAÇÕES: [ΔVdc, Δωe]

s = tf('s');
tau = L/R;  % Constante de tempo

%--- MODELO DESACOPLADO ---
fprintf('\nMODELO DESACOPLADO:\n');

% G(1,1): Id / Vd
G_des(1,1) = (1/L) / (s + R/L);
fprintf('G(1,1) = Id/Vd = %.2f / (%.4f*s + 1)\n', (1/L)*tau, tau);

% G(1,2): Id / Vq (desacoplado = 0)
G_des(1,2) = 0;
fprintf('G(1,2) = Id/Vq = 0 (desacoplado)\n');

% G(2,1): Iq / Vd (desacoplado = 0)
G_des(2,1) = 0;
fprintf('G(2,1) = Iq/Vd = 0 (desacoplado)\n');

% G(2,2): Iq / Vq
G_des(2,2) = (1/L) / (s + R/L);
fprintf('G(2,2) = Iq/Vq = %.2f / (%.4f*s + 1)\n', (1/L)*tau, tau);

%--- MODELO ACOPLADO ---
fprintf('\nMODELO ACOPLADO:\n');

% Sistema MIMO completo
C = eye(2);
D = zeros(2,2);
sys_acoplado = ss(F_acoplado, G_acoplado, C, D);
G_acop = tf(sys_acoplado);

fprintf('G(1,1) = Id/Vd (com acoplamento):\n');
G_acop(1,1)
fprintf('G(1,2) = Id/Vq (acoplamento cruzado):\n');
G_acop(1,2)
fprintf('G(2,1) = Iq/Vd (acoplamento cruzado):\n');
G_acop(2,1)
fprintf('G(2,2) = Iq/Vq (com acoplamento):\n');
G_acop(2,2)

%% MODELAGEM DE PERTURBAÇÕES
%==========================================================================
fprintf('\n=== PERTURBAÇÕES ===\n');

%--- Perturbação 1: Variação de Vdc ---
fprintf('\nPerturbação 1: ΔVdc (variação da tensão CC)\n');
fprintf('Efeito: Limita a máxima tensão dq aplicável\n');
fprintf('Modelo: Restrições do GPC devem considerar Vd_max = f(Vdc)\n');
fprintf('Não entra diretamente em Gq, mas afeta restrições\n');

%--- Perturbação 2: Variação de ωe ---
fprintf('\nPerturbação 2: Δωe (variação da frequência elétrica)\n');
fprintf('Efeito: Altera acoplamento cruzado entre eixos dq\n\n');

% No modelo acoplado, ωe aparece explicitamente em F
% Variação Δωe afeta os termos fora da diagonal

% Gq(1,1): efeito de Δωe em Id (via acoplamento com Iq)
% Linearizando: ΔId ≈ (∂Id/∂ωe)*Δωe
% Do modelo: dId/dt = ... + ωe*Iq
% Portanto: Δ(dId/dt) = Iq_ref * Δωe

Gq(1,1) = 0;  % ΔVdc não afeta diretamente Id
Gq(1,2) = x_ref(2) * (1/(s + R/L));  % Δωe → Id (via Iq)

fprintf('Gq(1,1) = ΔId/ΔVdc = 0 (efeito indireto via restrições)\n');
fprintf('Gq(1,2) = ΔId/Δωe ≈ %.4f / (%.4f*s + 1)\n', ...
        x_ref(2)*tau, tau);

% Gq(2,1): efeito de ΔVdc em Iq
Gq(2,1) = 0;

% Gq(2,2): efeito de Δωe em Iq (via acoplamento com Id)
Gq(2,2) = -x_ref(1) * (1/(s + R/L));

fprintf('Gq(2,1) = ΔIq/ΔVdc = 0 (efeito indireto via restrições)\n');
fprintf('Gq(2,2) = ΔIq/Δωe ≈ %.4f / (%.4f*s + 1)\n', ...
        -x_ref(1)*tau, tau);

%% DISCRETIZAÇÃO PARA GPC
%==========================================================================
fprintf('\n=== DISCRETIZAÇÃO ===\n');

% Escolher modelo (desacoplado ou acoplado)
USAR_ACOPLADO = true;

if USAR_ACOPLADO
    fprintf('Usando MODELO ACOPLADO para GPC\n');
    Gz = c2d(G_acop, Ts, 'zoh');
    A = A_acoplado;
    B = B_acoplado;
else
    fprintf('Usando MODELO DESACOPLADO para GPC\n');
    Gz = c2d(G_des, Ts, 'zoh');
    A = A_desacoplado;
    B = B_desacoplado;
end

Gzq = c2d(Gq, Ts, 'zoh');

fprintf('\nModelos discretizados (Ts = %.1f ms)\n', Ts*1000);

%% ANÁLISE DE GANHOS DC
%==========================================================================
fprintf('\n=== GANHOS DC (Estado Estacionário) ===\n');

Gz_dc = dcgain(Gz);
fprintf('Matriz de ganhos DC [Id,Iq]/[Vd,Vq]:\n');
disp(Gz_dc);

fprintf('\nInterpretação física:\n');
fprintf('Em regime permanente (DC):\n');
fprintf('  Id_ss = %.4f * Vd + %.4f * Vq [A/V]\n', Gz_dc(1,1), Gz_dc(1,2));
fprintf('  Iq_ss = %.4f * Vd + %.4f * Vq [A/V]\n', Gz_dc(2,1), Gz_dc(2,2));
fprintf('\nGanho esperado (R): 1/R = %.4f [A/V]\n', 1/R);

%% VALIDAÇÃO: RESPOSTA AO DEGRAU
%==========================================================================
figure('Name', 'Validação: Vd,Vq → Id,Iq', 'Position', [100 100 1200 600]);

t_sim = 0:Ts:0.05;  % 50 ms

subplot(1,2,1);
step(Gz(1,1), Gz(1,2), t_sim);
title('Resposta de I_d a degraus unitários em V_d e V_q');
xlabel('Tempo [s]'); ylabel('I_d [A]');
legend('V_d → I_d', 'V_q → I_d');
grid on;

subplot(1,2,2);
step(Gz(2,1), Gz(2,2), t_sim);
title('Resposta de I_q a degraus unitários em V_d e V_q');
xlabel('Tempo [s]'); ylabel('I_q [A]');
legend('V_d → I_q', 'V_q → I_q');
grid on;

%% ANÁLISE RGA (ACOPLAMENTO)
%==========================================================================
if USAR_ACOPLADO
    fprintf('\n=== ANÁLISE DE ACOPLAMENTO (RGA) ===\n');
    
    RGA = Gz_dc .* inv(Gz_dc)';
    fprintf('RGA (Relative Gain Array):\n');
    disp(RGA);
    
    fprintf('Diagonal RGA: [%.4f, %.4f]\n', RGA(1,1), RGA(2,2));
    
    if abs(RGA(1,1) - 1) < 0.1
        fprintf('→ Acoplamento FRACO: Controle SISO pode ser suficiente\n');
    else
        fprintf('→ Acoplamento FORTE: Controle MIMO é necessário\n');
    end
end

%% RESUMO FINAL
%==========================================================================
fprintf('\n========================================\n');
fprintf('MODELO GPC MIMO - CONFIGURAÇÃO CORRETA\n');
fprintf('========================================\n');
fprintf('VARIÁVEIS:\n');
fprintf('  • Saídas controladas (x): [Id, Iq] [A]\n');
fprintf('  • Entradas manipuladas (u): [Vd, Vq] [V]\n');
fprintf('  • Perturbações (z): [ΔVdc, Δωe] [V, rad/s]\n');
fprintf('\nCAMADA DE MODULAÇÃO PWM:\n');
fprintf('  • GPC calcula: Vd*, Vq* [V]\n');
fprintf('  • Transformações: dq → αβ → abc\n');
fprintf('  • PWM converte: Va*,Vb*,Vc* → da,db,dc\n');
fprintf('  • Relação: da = (Va* + Vdc/2) / Vdc\n');
fprintf('\nMODELO UTILIZADO:\n');
if USAR_ACOPLADO
    fprintf('  ✓ Acoplado (considera ωe explicitamente)\n');
else
    fprintf('  ✓ Desacoplado (ωe como perturbação)\n');
end
fprintf('  • Constante de tempo: τ = %.2f ms\n', tau*1000);
fprintf('  • Ganho DC: 1/R = %.4f A/V\n', 1/R);
fprintf('========================================\n');

%% SALVAR
save('modelo_gpc_correto.mat', 'Gz', 'Gzq', 'A', 'B', ...
     'Ts', 'R', 'L', 'Vdc', 'omega_e', 'x_ref', 'm', 'n', 'mq');
```

---

## 7. Resumo: Por Que Vd, Vq São as Entradas?

### **Argumento Físico**
```
Lei de Ohm + Lei de Faraday:
V = R*I + L*(dI/dt)

Resolvendo para dI/dt:
dI/dt = (V - R*I) / L

→ V é a causa (entrada)
→ I é o efeito (saída)
```

### **Argumento de Controle**
- **GPC decide**: "Que tensão aplicar para obter a corrente desejada?"
- **Não decide**: "Que duty cycle usar?" (isso é problema do PWM)

### **Hierarquia Completa**
```
GPC: Calcula Vd*, Vq* (referências de tensão)
  ↓
Transformações: dq→αβ→abc
  ↓
PWM: Converte Va*,Vb*,Vc* em pulsos
  ↓
Conversor: Aplica tensões ao circuito RL
  ↓
Sensores: Medem Id, Iq (realimentação)
