function [v_abc, duty, t_pwm] = applyPWM(in, out, i_init, i_sim, ...
                                           Ts_pwm, Ts_control, t_control, ...
                                           fs, Vdc, Nu, N1, N2, delta, lambda)
% =========================================================================
% applyPWM — SPWM Trifásico 2 Níveis
%
% Base teórica (conforme texto da dissertação):
%   u_x ∈ {0,1}  : posição da chave na fase x ∈ {a,b,c}
%   V_x = Vd/2 · u_x  : tensão de saída por fase (eq. 1 do texto)
%   8 estados de comutação → 7 vetores de tensão
%   Vetores nulos: (0,0,0)=V0 e (1,1,1)=V7
%
% SAÍDAS:
%   v_abc : struct — tensões geradas pelo inversor
%   duty  : struct — estados de chaveamento u_x ∈ {0,1}
%   t_pwm : vetor de tempo na resolução PWM [1 x n_pwm]
% =========================================================================

% -------------------------------------------------------------------------
% PASSO 1 — JANELA TEMPORAL E ÍNDICES
% -------------------------------------------------------------------------
% i_init: índice de início no vetor t_control
% i_sim:  índice atual (fim da janela desta chamada)
%
% O GPC chama applyPWM incrementalmente a cada Ts_control.
% A função opera apenas na janela [i_init, i_sim] — sem reprocessar
% instantes anteriores. Isso espelha a arquitetura real do DSP:
% a ISR de controle dispara applyPWM apenas para o período atual.

k_inicio = max(1, i_init);
k_fim    = min(i_sim, length(t_control));

t_ctrl_janela = t_control(k_inicio:k_fim);   % [1 x n_janela]
n_janela      = length(t_ctrl_janela);

% Razão de períodos: quantas amostras PWM por ciclo de controle
% Deve ser inteiro — sincronismo DSP/modulador
ratio = round(Ts_control / Ts_pwm);

% -------------------------------------------------------------------------
% PASSO 2 — EIXO DE TEMPO PWM
% -------------------------------------------------------------------------
% Resolução fina: passo Ts_pwm dentro da janela de controle.
% Cada intervalo [k, k+1] do controlador contém exatamente
% "ratio" amostras PWM — correspondendo a ratio ciclos de portadora.

t_pwm = t_ctrl_janela(1) : Ts_pwm : t_ctrl_janela(end);
n_pwm = length(t_pwm);

% -------------------------------------------------------------------------
% PASSO 3 — EXTRAÇÃO DAS REFERÊNCIAS DO GPC (in.u)
% -------------------------------------------------------------------------
% in.u [n_total x 3]: sinal de controle ótimo calculado pelo GPC
% Colunas: [v_a*  v_b*  v_c*] — referências de tensão por fase
%
% O GPC minimiza:
%   J = lambda·||y_pred - w||² + delta·||Δu||²
% e entrega v_abc* que applyPWM deve sintetizar via PWM.
%
% N1, N2: horizontes de predição (usados pelo GPC upstream)
% Nu:     horizonte de controle  (idem)
% Aqui apenas documentamos — applyPWM não recalcula o GPC.

if ~isfield(in,'u') || size(in.u,2) < 3
    error('applyPWM: in.u deve ser [n x 3] com colunas [va* vb* vc*]');
end

v_ref = in.u(k_inicio:k_fim, 1:3);   % [n_janela x 3]
%  v_ref(:,1) = v_a*  referência fase A
%  v_ref(:,2) = v_b*  referência fase B
%  v_ref(:,3) = v_c*  referência fase C

% -------------------------------------------------------------------------
% PASSO 4 — TABELA DOS 7 VETORES DE TENSÃO (2 Níveis)
% -------------------------------------------------------------------------
% Conforme texto: 8 estados → 7 vetores distintos
%
% Notação: (u_a, u_b, u_c) com u_x ∈ {0,1}
%
% Vetor | Estado     | v_ag         | v_bg         | v_cg
% ------|------------|--------------|--------------|-------------
%  V0   | (0, 0, 0)  |      0       |      0       |      0       ← nulo
%  V1   | (1, 0, 0)  |  +2Vdc/3    |  -Vdc/3      |  -Vdc/3
%  V2   | (1, 1, 0)  |  +Vdc/3     |  +Vdc/3      |  -2Vdc/3
%  V3   | (0, 1, 0)  |  -Vdc/3     |  +2Vdc/3     |  -Vdc/3
%  V4   | (0, 1, 1)  |  -2Vdc/3    |  +Vdc/3      |  +Vdc/3
%  V5   | (0, 0, 1)  |  -Vdc/3     |  -Vdc/3      |  +2Vdc/3
%  V6   | (1, 0, 1)  |  +Vdc/3     |  -2Vdc/3     |  +Vdc/3
%  V7   | (1, 1, 1)  |      0       |      0       |      0       ← nulo
%
% V0 e V7: vetores nulos — polo conectado ao mesmo barramento nas 3 fases
% Ponto na origem do plano αβ (conforme Figura 4 do texto)

% Tabela para log/debug [8 x 3]: cada linha = (u_a, u_b, u_c)
estados_comutacao = [0 0 0;   % V0 — nulo
                     1 0 0;   % V1
                     1 1 0;   % V2
                     0 1 0;   % V3
                     0 1 1;   % V4
                     0 0 1;   % V5
                     1 0 1;   % V6
                     1 1 1];  % V7 — nulo

% -------------------------------------------------------------------------
% PASSO 5 — PORTADORA TRIANGULAR SIMÉTRICA
% -------------------------------------------------------------------------
% A portadora triangular define os instantes de chaveamento.
% Amplitude ±1 (normalizada) — referência será dividida por Vdc/2.
%
% Equação da portadora triangular simétrica:
%   carrier(t) = (2/π)·arcsin(sin(2π·fs·t))   [forma analítica]
%   Implementação: sawtooth(2π·fs·t, 0.5)      [MATLAB — idêntico]
%
% Por que simétrica (width=0.5)?
%   - Pulsos de chaveamento centrados no período Ts_pwm
%   - Cancela harmônicos de ordem par (simetria de meia-onda)
%   - Harmônicos concentrados em fs, 2fs, 3fs... (fácil filtragem)
%   - Minimiza ripple de corrente para dado fs

carrier = sawtooth(2*pi*fs*t_pwm, 0.5);   % [1 x n_pwm], valores ∈ [-1,+1]

% -------------------------------------------------------------------------
% PASSO 6 — INTERPOLAÇÃO ZOH (Zero-Order Hold)
% -------------------------------------------------------------------------
% v_ref existe em t_ctrl_janela (passo Ts_control — mais lento)
% carrier existe em t_pwm       (passo Ts_pwm     — mais rápido)
%
% ZOH 'previous': mantém v_ref constante entre atualizações do GPC.
% Fisicamente: o registrador de comparação do timer PWM no DSP
% só é atualizado quando a ISR do controlador executa (Ts_control).
% Entre atualizações, o duty permanece fixo — comportamento ZOH exato.

va_pwm = interp1(t_ctrl_janela, v_ref(:,1), t_pwm, 'previous','extrap');
vb_pwm = interp1(t_ctrl_janela, v_ref(:,2), t_pwm, 'previous','extrap');
vc_pwm = interp1(t_ctrl_janela, v_ref(:,3), t_pwm, 'previous','extrap');
% Resultado: [1 x n_pwm] — cada valor repetido "ratio" vezes

% -------------------------------------------------------------------------
% PASSO 7 — TENSÃO DE SAÍDA CONFORME EQ. (1) DO TEXTO
% -------------------------------------------------------------------------
% O texto define:
%   V_x = Vd/2 · u_x    (equação 1)
%
% onde u_x ∈ {0,1} é a posição da chave (duty binário).
%
% Relação entre referência analógica e duty binário (via SPWM):
%   u_x = 1  quando  v_x_ref / (Vdc/2) >= carrier   → V_x = Vdc/2
%   u_x = 0  quando  v_x_ref / (Vdc/2) <  carrier   → V_x = 0
%
% Portanto:
%   Normalização: Ma_x = v_x_ref / (Vdc/2)  ∈ [-1, +1]
%   Comparação:   u_x  = (Ma_x >= carrier)  ∈ {0, 1}
%   Tensão polo:  V_x  = (Vdc/2) · u_x     ∈ {0, Vdc/2}
%
% ATENÇÃO — diferença de referencial:
%   Eq.(1) do texto: V_x referida ao NEUTRO do barramento (±Vdc/2)
%   Implementação:   v_XN referida ao NEGATIVO do barramento {0, Vdc}
%   Conversão:       v_Xg (fase-neutro carga) remove tensão de modo comum

% Normalização — índice de modulação por fase
Ma_a = va_pwm / (Vdc/2);   % [1 x n_pwm]
Ma_b = vb_pwm / (Vdc/2);
Ma_c = vc_pwm / (Vdc/2);

% Aviso de sobremodulação
Ma_max = max(abs([Ma_a, Ma_b, Ma_c]));
if Ma_max > 1.0
    warning('applyPWM: Sobremodulacao | Ma_max=%.3f > 1.0 | GPC saturou referencia', Ma_max);
end

% -------------------------------------------------------------------------
% PASSO 8 — GERAÇÃO DE u_x ∈ {0,1} — ESTADOS DE CHAVEAMENTO
% -------------------------------------------------------------------------
% Implementação direta da lógica SPWM:
%   u_x = 1  ↔  S_superior FECHADA  ↔  fase conectada a +Vdc/2
%   u_x = 0  ↔  S_inferior FECHADA  ↔  fase conectada a -Vdc/2
%
% Comparação: referência normalizada vs portadora triangular
% Resultado: u_x ∈ {0,1} conforme notação do texto

u_a = double(Ma_a >= carrier);   % [1 x n_pwm]  u_a ∈ {0,1}
u_b = double(Ma_b >= carrier);   % [1 x n_pwm]  u_b ∈ {0,1}
u_c = double(Ma_c >= carrier);   % [1 x n_pwm]  u_c ∈ {0,1}

% -------------------------------------------------------------------------
% PASSO 9 — TENSÃO DE POLO: V_x = Vdc/2 · u_x  (Equação 1 do texto)
% -------------------------------------------------------------------------
% Aplicação direta da equação (1):
%   V_a = Vdc/2 · u_a
%   V_b = Vdc/2 · u_b
%   V_c = Vdc/2 · u_c
%
% Resultado: tensão referida ao ponto médio do barramento CC
%   u_x = 1 → V_x = +Vdc/2  (polo em +)
%   u_x = 0 → V_x =  0      (polo em -, referido ao neutro = -Vdc/2)
%
% Para cálculo das tensões de fase usaremos v_XN = u_x · Vdc
% (referida ao negativo) — algebricamente equivalente após subtração
% da tensão de modo comum

V_a = (Vdc/2) * u_a;   % equação (1) — [1 x n_pwm] ∈ {0, Vdc/2}
V_b = (Vdc/2) * u_b;
V_c = (Vdc/2) * u_c;

% Versão referida ao negativo (para cálculo de v_fase):
% v_XN = u_x · Vdc  (desloca referencial de -Vdc/2 para 0)
v_AN = u_a * Vdc;   % [1 x n_pwm] ∈ {0, Vdc}
v_BN = u_b * Vdc;
v_CN = u_c * Vdc;

% -------------------------------------------------------------------------
% PASSO 10 — IDENTIFICAÇÃO DO VETOR DE TENSÃO ATIVO
% -------------------------------------------------------------------------
% A cada instante PWM, o estado (u_a, u_b, u_c) define um dos
% 8 estados → 7 vetores de tensão distintos.
%
% Vetores nulos (V0, V7): u_a=u_b=u_c → todos os polos no mesmo barramento
%   V0=(0,0,0): todas fases em -Vdc/2 → tensão diferencial = 0
%   V7=(1,1,1): todas fases em +Vdc/2 → tensão diferencial = 0
%
% Identificação numérica: índice = u_a·4 + u_b·2 + u_c·1  (binário)
%   (0,0,0)=0→V0  (1,0,0)=4→V1  (1,1,0)=6→V2  (0,1,0)=2→V3
%   (0,1,1)=3→V4  (0,0,1)=1→V5  (1,0,1)=5→V6  (1,1,1)=7→V7

vetor_idx = u_a*4 + u_b*2 + u_c*1;   % [1 x n_pwm] ∈ {0,1,2,3,4,5,6,7}

% Contagem de tempo em cada vetor (para análise de dissertação)
tempo_vetor = zeros(1,8);
for v = 0:7
    tempo_vetor(v+1) = sum(vetor_idx == v) * Ts_pwm;
end

% -------------------------------------------------------------------------
% PASSO 11 — TENSÕES DE FASE (phase-to-neutral)
% -------------------------------------------------------------------------
% Carga trifásica equilibrada com neutro flutuante:
%
%   Tensão de modo comum: v_gN = (v_AN + v_BN + v_CN) / 3
%
%   v_ag = v_AN - v_gN = (2·v_AN - v_BN - v_CN) / 3
%   v_bg = v_BN - v_gN = (2·v_BN - v_AN - v_CN) / 3
%   v_cg = v_CN - v_gN = (2·v_CN - v_AN - v_BN) / 3
%
% Equivalentemente usando V_x da equação (1):
%   v_ag = V_a - (V_a + V_b + V_c)/3  [referido ao ponto médio]
%
% Níveis discretos de v_ag no conversor 2 níveis:
%   {-2Vdc/3, -Vdc/3, 0, +Vdc/3, +2Vdc/3}

v_ag = (2*v_AN - v_BN - v_CN) / 3;
v_bg = (2*v_BN - v_AN - v_CN) / 3;
v_cg = (2*v_CN - v_AN - v_BN) / 3;

% Verificação: v_ag + v_bg + v_cg = 0 (equilíbrio trifásico)
assert(max(abs(v_ag + v_bg + v_cg)) < 1e-9, ...
       'applyPWM: ERRO — tensoes de fase nao somam zero');

% -------------------------------------------------------------------------
% PASSO 12 — TENSÕES DE LINHA (line-to-line)
% -------------------------------------------------------------------------
% v_ab = v_AG - v_BG = v_AN - v_BN
%
% 2 Níveis: v_linha ∈ {-Vdc, 0, +Vdc}  (3 níveis de tensão de linha)
% Pico teórico SPWM (Ma=1): v_ab_pico = Vdc·√3/2 ≈ 0.866·Vdc

v_ab = v_AN - v_BN;
v_bc = v_BN - v_CN;
v_ca = v_CN - v_AN;

% -------------------------------------------------------------------------
% PASSO 13 — DUTY MÉDIO POR PERÍODO DE CONTROLE
% -------------------------------------------------------------------------
% O modelo médio do conversor (usado pelo GPC como planta interna):
%   <v_ag> ≈ (d_a - 0.5)·Vdc  onde d_a = duty médio de u_a
%
% Útil para verificar coerência entre referência GPC e saída PWM.
% Discrepância indica saturação ou não-linearidade da modulação.

duty_media_a = zeros(1, n_janela);
duty_media_b = zeros(1, n_janela);
duty_media_c = zeros(1, n_janela);

for k = 1:n_janela
    idx_bloco = (k-1)*ratio + 1 : min(k*ratio, n_pwm);
    duty_media_a(k) = mean(u_a(idx_bloco));
    duty_media_b(k) = mean(u_b(idx_bloco));
    duty_media_c(k) = mean(u_c(idx_bloco));
end

% -------------------------------------------------------------------------
% PASSO 14 — MONTAGEM DAS STRUCTS DE SAÍDA
% -------------------------------------------------------------------------

% Tensões do inversor
v_abc.fase   = [v_ag; v_bg; v_cg];   % [3 x n_pwm]  fase-neutro
v_abc.linha  = [v_ab; v_bc; v_ca];   % [3 x n_pwm]  linha-linha
v_abc.polo   = [V_a;  V_b;  V_c ];   % [3 x n_pwm]  eq.(1): Vdc/2·u_x
v_abc.modo_comum = (v_AN + v_BN + v_CN)/3; % [1 x n_pwm]

% Estados de chaveamento
duty.u       = [u_a; u_b; u_c];         % [3 x n_pwm]  u_x ∈ {0,1}
duty.media   = [duty_media_a; ...
                duty_media_b; ...
                duty_media_c];           % [3 x n_janela]
duty.Ma      = [Ma_a; Ma_b; Ma_c];      % [3 x n_pwm]
duty.vetor   = vetor_idx;               % [1 x n_pwm]  ∈ {0..7}
duty.t_vetor = tempo_vetor;             % [1 x 8]  tempo em cada vetor (s)

% =========================================================================
% PASSO 15 — PLOTS INTERNOS
% =========================================================================

t_ms      = t_pwm * 1e3;
t_ctrl_ms = t_ctrl_janela * 1e3;

% Janela de zoom: últimos 6 períodos PWM (3 completos visíveis)
n_zoom = min(6*ratio, n_pwm);
idx_z  = (n_pwm - n_zoom + 1) : n_pwm;
t_us   = t_pwm(idx_z) * 1e6;

% Janela de regime: últimos 2 ciclos de rede (40 ms a 50 Hz)
t_reg  = max(t_ms(1), t_ms(end) - 2*20);
idx_r  = t_ms >= t_reg;
idx_rc = t_ctrl_ms >= t_reg;

%-- Plot 1: Portadora + Referências + u_x ---------------------------------
figure('Name','[PWM 2-Nível] Portadora, Referências e Chaveamento', ...
       'Color','w','Position',[40 560 1000 480]);

subplot(3,1,1)
plot(t_us, carrier(idx_z), 'k-', 'LineWidth',2.0); hold on
plot(t_us, Ma_a(idx_z),    'b-', 'LineWidth',1.6)
plot(t_us, Ma_b(idx_z),    'r-', 'LineWidth',1.6)
plot(t_us, Ma_c(idx_z),    'g-', 'LineWidth',1.6)
yline( 1,'k:'); yline(-1,'k:')
legend('Portadora','Ma_a','Ma_b','Ma_c','Location','best','FontSize',8)
ylabel('Amplitude Norm.'); grid on; ylim([-1.3 1.3])
title(sprintf('SPWM 2 Níveis | f_s=%gkHz | V_{dc}=%gV | Ma_{max}=%.3f', ...
              fs/1e3, Vdc, Ma_max))

subplot(3,1,2)
% Estados u_x empilhados para visualização
stairs(t_us, u_a(idx_z),   'b-', 'LineWidth',1.5); hold on
stairs(t_us, u_b(idx_z)+2, 'r-', 'LineWidth',1.5)
stairs(t_us, u_c(idx_z)+4, 'g-', 'LineWidth',1.5)
yticks([0 1 2 3 4 5]); yticklabels({'0','1','0','1','0','1'})
legend('u_a','u_b','u_c','Location','best','FontSize',8)
ylabel('u_x \in {0,1}'); grid on
title('Estados de Chaveamento u_x \in {0,1} por Fase (empilhados)')

subplot(3,1,3)
% Tensão de polo: V_x = Vdc/2 · u_x  (equação 1 do texto)
stairs(t_us, V_a(idx_z),          'b-', 'LineWidth',1.4); hold on
stairs(t_us, V_b(idx_z)+Vdc*0.6,  'r-', 'LineWidth',1.4)
stairs(t_us, V_c(idx_z)+Vdc*1.2,  'g-', 'LineWidth',1.4)
ylabel('V_x = V_{dc}/2·u_x (V)'); xlabel('Tempo (µs)'); grid on
title('Tensão de Polo por Fase — Eq.(1): V_x = V_d/2 · u_x')
legend('V_a','V_b','V_c','Location','best','FontSize',8)

%-- Plot 2: Vetores de Tensão Ativos (sequência temporal) -----------------
figure('Name','[PWM 2-Nível] Sequência de Vetores','Color','w', ...
       'Position',[40 60 1000 300]);

nomes_vetores = {'V0(000)','V5(001)','V3(010)','V4(011)', ...
                 'V1(100)','V6(101)','V2(110)','V7(111)'};
cores_v = [0.5 0.5 0.5; 0.2 0.6 1; 0.2 0.8 0.2; 0.1 0.5 0.1;
           1 0.3 0.3;   0.8 0.1 0.8; 1 0.6 0; 0.3 0.3 0.3];

stairs(t_ms(idx_z), vetor_idx(idx_z), 'k-', 'LineWidth', 1.5); hold on
for v = 0:7
    mask = vetor_idx(idx_z) == v;
    t_v  = t_us(mask);
    if any(mask)
        scatter(t_v, v*ones(size(t_v)), 20, cores_v(v+1,:), 'filled')
    end
end
yticks(0:7)
yticklabels({'V0(000)','V1(100) ← texto','V2(110)','V3(010)', ...
             'V4(011)','V5(001)','V6(101)','V7(111)'})
xlabel('Tempo (µs)'); ylabel('Vetor Ativo')
title('Sequência de Vetores de Tensão — 8 Estados de Comutação (2 Níveis)')
grid on

%-- Plot 3: Tensões de Fase e Polo ----------------------------------------
figure('Name','[PWM 2-Nível] Tensões de Fase','Color','w', ...
       'Position',[60 60 1000 600]);

subplot(2,1,1)
plot(t_ms(idx_r), v_ag(idx_r), 'b', 'LineWidth',0.8); hold on
plot(t_ms(idx_r), v_bg(idx_r), 'r', 'LineWidth',0.8)
plot(t_ms(idx_r), v_cg(idx_r), 'g', 'LineWidth',0.8)
plot(t_ctrl_ms(idx_rc), v_ref(idx_rc,1), 'b--', 'LineWidth',1.8)
plot(t_ctrl_ms(idx_rc), v_ref(idx_rc,2), 'r--', 'LineWidth',1.8)
plot(t_ctrl_ms(idx_rc), v_ref(idx_rc,3), 'g--', 'LineWidth',1.8)
yline( 2*Vdc/3,'k:','+2V_{dc}/3','LabelHorizontalAlignment','left')
yline(-2*Vdc/3,'k:','-2V_{dc}/3','LabelHorizontalAlignment','left')
ylabel('Tensão (V)'); grid on
legend('v_{ag}','v_{bg}','v_{cg}','v_a^*','v_b^*','v_c^*', ...
       'Location','best','FontSize',7)
title(sprintf('Tensões de Fase | Níveis: {0, ±V_{dc}/3, ±2V_{dc}/3} | V_{dc}=%gV',Vdc))
xlim([t_reg t_ms(end)])

subplot(2,1,2)
plot(t_ms(idx_r), v_ab(idx_r), 'b', 'LineWidth',0.8); hold on
plot(t_ms(idx_r), v_bc(idx_r), 'r', 'LineWidth',0.8)
plot(t_ms(idx_r), v_ca(idx_r), 'g', 'LineWidth',0.8)
yline( Vdc,'k:','+V_{dc}','LabelHorizontalAlignment','left')
yline(-Vdc,'k:','-V_{dc}','LabelHorizontalAlignment','left')
ylabel('Tensão (V)'); xlabel('Tempo (ms)'); grid on
legend('v_{ab}','v_{bc}','v_{ca}','Location','best')
title('Tensões de Linha | 3 Níveis: {-V_{dc}, 0, +V_{dc}}')
xlim([t_reg t_ms(end)])

%-- Plot 4: Tempo em cada vetor (pizza + barra) ---------------------------
figure('Name','[PWM 2-Nível] Distribuição dos Vetores','Color','w', ...
       'Position',[500 400 900 380]);

subplot(1,2,1)
nomes_v = {'V0(000)','V1(100)','V2(110)','V3(010)', ...
           'V4(011)','V5(001)','V6(101)','V7(111)'};
t_v_ms = tempo_vetor * 1e3;
idx_ativos = t_v_ms > 0;
pie(t_v_ms(idx_ativos))
legend(nomes_v(idx_ativos), 'Location','bestoutside','FontSize',7)
title('Distribuição Temporal dos Vetores')

subplot(1,2,2)
bar(0:7, t_v_ms, 'FaceColor',[0.2 0.5 0.8])
xticks(0:7); xticklabels(nomes_v)
xtickangle(30)
xlabel('Vetor'); ylabel('Tempo (ms)')
title('Tempo em cada Vetor de Tensão')
yline(sum(t_v_ms)/8,'r--','Média','LabelHorizontalAlignment','right')
grid on

% Marca vetores nulos (V0 e V7) conforme texto
hold on
bar([0 7], t_v_ms([1 8]), 'FaceColor',[0.8 0.2 0.2])
legend('Vetores ativos','Vetores nulos (V0,V7)','Location','best','FontSize',8)

% -------------------------------------------------------------------------
% RELATÓRIO CONSOLE
% -------------------------------------------------------------------------
fprintf('\n=== applyPWM — Conversor 2 Níveis ===\n');
fprintf('  Equação (1): V_x = Vdc/2 · u_x\n');
fprintf('  u_x ∈ {0,1},  x ∈ {a,b,c}\n');
fprintf('  Janela     : k=%d → k=%d  (%d instantes ctrl)\n', k_inicio, k_fim, n_janela);
fprintf('  n_pwm      : %d amostras  (ratio=%d por ctrl)\n', n_pwm, ratio);
fprintf('  Ma_max     : %.4f', Ma_max);
if Ma_max <= 1.0
    fprintf('  [LINEAR OK]\n');
else
    fprintf('  [SOBREMODULACAO]\n');
end
fprintf('  Vetores    : V0=%dms  V7=%dms  (nulos)\n', ...
        round(tempo_vetor(1)*1e3), round(tempo_vetor(8)*1e3));
fprintf('  v_fase max : ±%.1f V  (= ±2Vdc/3)\n', 2*Vdc/3);
fprintf('  v_linha max: ±%.1f V  (= ±Vdc)\n',    Vdc);
fprintf('=====================================\n\n');

end
