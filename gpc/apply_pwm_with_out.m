function [v_abc, duty, t_pwm] = applyPWM(in, out, i_init, i_sim, ...
                                           Ts_pwm, Ts_control, t_control, ...
                                           fs, Vdc, Nu, N1, N2, delta, lambda)
% =========================================================================
% applyPWM — Modulação SPWM para Conversor Trifásico 2 Níveis
%
% SAÍDAS:
%   v_abc  : struct com tensões geradas pelo inversor
%              v_abc.fase  [3 x n_pwm]  tensões de fase  (v_ag, v_bg, v_cg)
%              v_abc.linha [3 x n_pwm]  tensões de linha (v_ab, v_bc, v_ca)
%              v_abc.polo  [3 x n_pwm]  tensões de polo  (v_AN, v_BN, v_CN)
%
%   duty   : struct com sinais de chaveamento
%              duty.sinal  [3 x n_pwm]  duty cycle binário {0,1}
%              duty.media  [3 x n_ctrl] duty médio por período de controle
%              duty.Ma     [3 x n_pwm]  índice de modulação normalizado
%
%   t_pwm  : vetor de tempo na resolução PWM [1 x n_pwm] (s)
% =========================================================================

% -------------------------------------------------------------------------
% PASSO 1 — ÍNDICES E JANELA TEMPORAL
% -------------------------------------------------------------------------
% i_init: iteração inicial — define onde começa esta chamada do applyPWM
% i_sim:  iteração atual   — define onde termina
%
% Justificativa: o GPC chama applyPWM a cada Ts_control.
% A função opera apenas sobre a janela [i_init, i_sim] do vetor t_control.
% Isso permite que applyPWM seja chamada incrementalmente dentro do loop GPC
% sem reprocessar instantes anteriores.

k_inicio = max(1,     i_init);
k_fim    = min(i_sim, length(t_control));

% Instantes de controle desta janela
t_ctrl_janela = t_control(k_inicio : k_fim);   % [1 x n_janela]
n_janela      = length(t_ctrl_janela);

% -------------------------------------------------------------------------
% PASSO 2 — EIXO DE TEMPO PWM
% -------------------------------------------------------------------------
% O eixo PWM tem passo Ts_pwm (muito menor que Ts_control).
% Cada intervalo de controle contém exatamente ratio amostras PWM.
%
% ratio = Ts_control / Ts_pwm  (deve ser inteiro — verificado no Passo 1)
%
% Exemplo: Ts_control = 200µs, Ts_pwm = 100µs → ratio = 2
%   t_control: [0,  200µs, 400µs, ...]
%   t_pwm:     [0, 100µs, 200µs, 300µs, 400µs, ...]

t_inicio = t_ctrl_janela(1);
t_fim    = t_ctrl_janela(end);

t_pwm = t_inicio : Ts_pwm : t_fim;   % [1 x n_pwm]
n_pwm = length(t_pwm);

% Verificação do ratio
ratio = round(Ts_control / Ts_pwm);

% -------------------------------------------------------------------------
% PASSO 3 — EXTRAÇÃO DAS REFERÊNCIAS DE TENSÃO (in.u)
% -------------------------------------------------------------------------
% in.u contém o sinal de controle ótimo calculado pelo GPC.
% Formato esperado: in.u [n_total_ctrl x 3] → colunas: [v_a*, v_b*, v_c*]
%
% Extraímos apenas a janela correspondente a [k_inicio, k_fim].
%
% Justificativa: o GPC pode calcular horizontes futuros (N2 passos),
% mas applyPWM aplica apenas os instantes já decorridos.

if ~isfield(in, 'u')
    error('applyPWM: in.u nao encontrado. GPC deve preencher in.u [n x 3].');
end

if size(in.u, 2) < 3
    error('applyPWM: in.u deve ter 3 colunas [v_a* v_b* v_c*].');
end

% Referências de tensão na janela atual [n_janela x 3]
v_ref = in.u(k_inicio:k_fim, 1:3);
%        col1=v_a*  col2=v_b*  col3=v_c*

% -------------------------------------------------------------------------
% PASSO 4 — PORTADORA TRIANGULAR SIMÉTRICA
% -------------------------------------------------------------------------
% A portadora triangular é o elemento central do SPWM.
%
% CONVERSOR 2 NÍVEIS — características da portadora:
%   - Frequência : fs (igual à frequência de chaveamento)
%   - Amplitude  : ±1 (normalizada)
%   - Simetria   : dupla (pico em 0.5 do período) → cancela harmônicos pares
%
% sawtooth(x, width):
%   width = 0.5 → triângulo simétrico (subida = descida = Ts_pwm/2)
%   width = 1   → dente de serra (apenas descida) — NÃO usar para SPWM
%
% Por que triangular simétrica?
%   No conversor 2 níveis, a portadora simétrica garante que os pulsos
%   de chaveamento sejam centrados no período → menor ripple de corrente
%   para a mesma frequência de chaveamento (equivalente ao SVPWM em
%   termos de conteúdo harmônico médio).

carrier = sawtooth(2*pi*fs*t_pwm, 0.5);
% Resultado: valores em [-1, +1], período = 1/fs = Ts_pwm

% -------------------------------------------------------------------------
% PASSO 5 — INTERPOLAÇÃO ZOH (Zero-Order Hold)
% -------------------------------------------------------------------------
% Problema: v_ref existe em t_ctrl_janela (passo Ts_control)
%           carrier existe em t_pwm       (passo Ts_pwm)
%           → precisamos de v_ref no mesmo passo de t_pwm
%
% Solução: ZOH com 'previous' — mantém o valor constante até
%          a próxima atualização do controlador.
%
% Por que 'previous' e não 'linear'?
%   No DSP real, o registrador de comparação do timer PWM é carregado
%   uma vez por período de controle e mantido fixo até a próxima
%   interrupção (ISR). Isso é exatamente um ZOH.
%   Interpolação linear implicaria variação contínua do duty — fisicamente
%   impossível sem hardware adicional.
%
% 'extrap': repete o último valor fora do intervalo (borda final)

va_pwm = interp1(t_ctrl_janela, v_ref(:,1), t_pwm, 'previous', 'extrap');
vb_pwm = interp1(t_ctrl_janela, v_ref(:,2), t_pwm, 'previous', 'extrap');
vc_pwm = interp1(t_ctrl_janela, v_ref(:,3), t_pwm, 'previous', 'extrap');
% Resultado: [1 x n_pwm] — cada valor repetido por "ratio" amostras

% -------------------------------------------------------------------------
% PASSO 6 — NORMALIZAÇÃO: ÍNDICE DE MODULAÇÃO
% -------------------------------------------------------------------------
% O comparador PWM opera com sinais normalizados em [-1, +1].
% Normalização: Ma_x = v_x* / (Vdc/2)
%
% CONVERSOR 2 NÍVEIS — limites de operação:
%   Região linear:      |Ma| ≤ 1.0   → tensão fundamental proporcional
%   Sobremodulação:     1.0 < |Ma| ≤ 3.24  → distorção harmônica crescente
%   Onda quadrada:      |Ma| → ∞     → máxima tensão fundamental = 4/π·Vdc/2
%
% Para SPWM linear trifásico:
%   Tensão de fase máxima (pico) = Vdc/2
%   Tensão de linha máxima (pico) = Vdc·√3/2 ≈ 0.866·Vdc
%
% Importante: o GPC deve limitar suas referências a ±Vdc/2
%   Se |v_ref| > Vdc/2 → Ma > 1 → sobremodulação

Ma_a = va_pwm / (Vdc/2);    % [1 x n_pwm]
Ma_b = vb_pwm / (Vdc/2);
Ma_c = vc_pwm / (Vdc/2);

% Verificação de sobremodulação
Ma_max = max(abs([Ma_a, Ma_b, Ma_c]));
if Ma_max > 1.0
    warning('applyPWM: Sobremodulacao detectada. Ma_max = %.3f > 1.0', Ma_max);
    warning('          Tensao de fase limitada a Vdc/2 = %.1f V', Vdc/2);
end

% -------------------------------------------------------------------------
% PASSO 7 — COMPARAÇÃO: GERAÇÃO DOS DUTY CYCLES
% -------------------------------------------------------------------------
% Lógica de chaveamento do inversor 2 níveis (por fase):
%
%   Ma_x(k) >= carrier(k)  →  duty_x = 1  →  S_superior FECHADA
%                                              S_inferior ABERTA
%                                              v_xN = Vdc
%
%   Ma_x(k) <  carrier(k)  →  duty_x = 0  →  S_superior ABERTA
%                                              S_inferior FECHADA
%                                              v_xN = 0
%
% Topologia 2 níveis:
%   +Vdc/2 ─── S_sup ─┐
%                      ├── X (saída de fase)
%   -Vdc/2 ─── S_inf ─┘
%
% Tensão de polo: v_XN ∈ {0, Vdc}  (referida ao negativo do barramento)
% Tensão de fase: v_Xg ∈ {-Vdc/3, -Vdc/6, +Vdc/6, +Vdc/3, +Vdc/2, -Vdc/2}
%                 (valores discretos dependem do estado das 3 fases)

duty_a = double(Ma_a >= carrier);   % [1 x n_pwm], valores em {0, 1}
duty_b = double(Ma_b >= carrier);
duty_c = double(Ma_c >= carrier);

% -------------------------------------------------------------------------
% PASSO 8 — TENSÕES DE POLO (Phase-to-DC-Negative)
% -------------------------------------------------------------------------
% v_XN: tensão entre o polo de fase X e o negativo do barramento CC.
%
% Conversor 2 níveis:
%   duty = 1 → chave superior fecha → polo conectado a +Vdc → v_XN = Vdc
%   duty = 0 → chave inferior fecha → polo conectado a GND  → v_XN = 0
%
% Nota: referido ao negativo do barramento, não ao neutro da carga.

v_AN = duty_a * Vdc;    % [1 x n_pwm], valores em {0, Vdc}
v_BN = duty_b * Vdc;
v_CN = duty_c * Vdc;

% -------------------------------------------------------------------------
% PASSO 9 — TENSÕES DE FASE (Phase-to-Neutral)
% -------------------------------------------------------------------------
% Para carga trifásica equilibrada com neutro flutuante:
%
%   v_Xg = v_XN - v_gN   (onde g = neutro da carga, N = negativo CC)
%
% Neutro flutuante: v_gN = (v_AN + v_BN + v_CN) / 3
%
% Substituindo:
%   v_ag = v_AN - (v_AN + v_BN + v_CN)/3 = (2·v_AN - v_BN - v_CN) / 3
%   v_bg = (2·v_BN - v_AN - v_CN) / 3
%   v_cg = (2·v_CN - v_AN - v_BN) / 3
%
% Verificação: v_ag + v_bg + v_cg = 0  (sistema equilibrado)
%
% Níveis de tensão de fase no conversor 2 níveis:
%   Estados possíveis de (duty_a, duty_b, duty_c):
%   (1,0,0) → v_ag = +2/3·Vdc   (1,1,0) → v_ag = +1/3·Vdc
%   (0,1,1) → v_ag = -2/3·Vdc   (0,0,1) → v_ag = -1/3·Vdc
%   (1,1,1) → v_ag = 0           (0,0,0) → v_ag = 0

v_ag = (2*v_AN - v_BN - v_CN) / 3;    % [1 x n_pwm]
v_bg = (2*v_BN - v_AN - v_CN) / 3;
v_cg = (2*v_CN - v_AN - v_BN) / 3;

% Verificação numérica do equilíbrio (tolerância por arredondamento float)
erro_equil = max(abs(v_ag + v_bg + v_cg));
if erro_equil > 1e-9
    warning('applyPWM: Desequilíbrio numérico nas tensões de fase: %.2e V', erro_equil);
end

% -------------------------------------------------------------------------
% PASSO 10 — TENSÕES DE LINHA (Line-to-Line)
% -------------------------------------------------------------------------
% v_ab = v_ag - v_bg = v_AN - v_BN
%
% Conversor 2 níveis — níveis de tensão de linha:
%   v_ab ∈ {-Vdc, 0, +Vdc}   (3 níveis — característica do 2 níveis)
%
% Pico teórico SPWM linear: v_ab_pico = Ma · Vdc · √3/2
% Para Ma = 1: v_ab_pico = 0.866 · Vdc

v_ab = v_AN - v_BN;    % equivalente a v_ag - v_bg
v_bc = v_BN - v_CN;
v_ca = v_CN - v_AN;

% -------------------------------------------------------------------------
% PASSO 11 — DUTY CYCLE MÉDIO POR PERÍODO DE CONTROLE
% -------------------------------------------------------------------------
% O duty médio é útil para:
%   (a) Verificar se o GPC está saturando a modulação
%   (b) Alimentar modelos de estado médio (averaged model)
%   (c) Comparar com a referência do MPC (deve ser ≈ Ma/2 + 0.5)
%
% duty_medio_x(k) = média de duty_x nos "ratio" amostras PWM
%                   dentro do k-ésimo período de controle

duty_media_a = zeros(1, n_janela);
duty_media_b = zeros(1, n_janela);
duty_media_c = zeros(1, n_janela);

for k = 1:n_janela
    idx_bloco = (k-1)*ratio + 1 : min(k*ratio, n_pwm);
    duty_media_a(k) = mean(duty_a(idx_bloco));
    duty_media_b(k) = mean(duty_b(idx_bloco));
    duty_media_c(k) = mean(duty_c(idx_bloco));
end

% -------------------------------------------------------------------------
% PASSO 12 — MONTAGEM DAS STRUCTS DE SAÍDA
% -------------------------------------------------------------------------

% Tensões geradas pelo inversor
v_abc.fase  = [v_ag;  v_bg;  v_cg ];   % [3 x n_pwm] tensões de fase
v_abc.linha = [v_ab;  v_bc;  v_ca ];   % [3 x n_pwm] tensões de linha
v_abc.polo  = [v_AN;  v_BN;  v_CN ];   % [3 x n_pwm] tensões de polo

% Sinais de chaveamento
duty.sinal  = [duty_a; duty_b; duty_c];            % [3 x n_pwm] binário
duty.media  = [duty_media_a; duty_media_b; duty_media_c]; % [3 x n_janela]
duty.Ma     = [Ma_a;  Ma_b;  Ma_c  ];              % [3 x n_pwm] normalizado

% -------------------------------------------------------------------------
% PASSO 13 — PLOTS INTERNOS
% -------------------------------------------------------------------------

t_ms    = t_pwm * 1e3;
t_ctrl_ms = t_ctrl_janela * 1e3;

% Zoom: últimos 5 períodos PWM (visualização da comutação)
n_zoom  = min(5*ratio, n_pwm);
idx_z   = (n_pwm - n_zoom + 1) : n_pwm;
t_us_z  = t_pwm(idx_z) * 1e6;   % µs

%-- Plot 1: Portadora vs Referências Normalizadas -------------------------
figure('Name','[PWM] Portadora vs Referências','Color','w', ...
       'Position',[50 550 950 380]);

subplot(2,1,1)
plot(t_us_z, carrier(idx_z),  'k-',  'LineWidth', 1.8); hold on
plot(t_us_z, Ma_a(idx_z),     'b-',  'LineWidth', 1.5)
plot(t_us_z, Ma_b(idx_z),     'r-',  'LineWidth', 1.5)
plot(t_us_z, Ma_c(idx_z),     'g-',  'LineWidth', 1.5)
yline( 1,'k:','LineWidth',0.8); yline(-1,'k:','LineWidth',0.8)
ylim([-1.3 1.3])
legend('Portadora',  'M_a = v_a^*/(V_{dc}/2)', 'M_b', 'M_c', ...
       'Location','best','FontSize',8)
title(sprintf('SPWM 2 Níveis — Portadora Triangular vs Referências | f_s = %g kHz', fs/1e3))
ylabel('Amplitude Normalizada'); grid on

subplot(2,1,2)
stairs(t_us_z, duty_a(idx_z),   'b-', 'LineWidth',1.4); hold on
stairs(t_us_z, duty_b(idx_z)+2, 'r-', 'LineWidth',1.4)
stairs(t_us_z, duty_c(idx_z)+4, 'g-', 'LineWidth',1.4)
yticks([0 1 2 3 4 5])
yticklabels({'0','1','0','1','0','1'})
ylabel('Duty  [a / b / c]'); xlabel('Tempo (µs)')
title('Duty Cycles — Chaveamento Binário {0,1} por Fase (empilhados)')
grid on

%-- Plot 2: Tensões de Polo -----------------------------------------------
figure('Name','[PWM] Tensões de Polo — v_XN','Color','w', ...
       'Position',[50 100 950 450]);

fases_polo  = {v_AN, v_BN, v_CN};
nomes_polo  = {'v_{AN}', 'v_{BN}', 'v_{CN}'};
cores_polo  = {'b','r','g'};

for f = 1:3
    subplot(3,1,f)
    stairs(t_ms(idx_z), fases_polo{f}(idx_z), ...
           cores_polo{f}, 'LineWidth', 1.3)
    ylabel(sprintf('%s (V)', nomes_polo{f}))
    ylim([-0.1*Vdc, 1.1*Vdc])
    yline(Vdc, 'k--', sprintf('V_{dc}=%gV',Vdc), 'LabelHorizontalAlignment','left')
    yline(0,   'k--')
    grid on
    if f==1
        title(sprintf('Tensões de Polo — Conversor 2 Níveis | V_{dc} = %g V', Vdc))
    end
    if f==3, xlabel('Tempo (ms)'); end
end

%-- Plot 3: Tensões de Fase e Linha ---------------------------------------
figure('Name','[PWM] Tensões de Fase e Linha','Color','w', ...
       'Position',[100 100 1000 600]);

% Janela maior: últimos 2 ciclos de rede
f_rede   = 50;
t_regime = max(t_ms(1), t_ms(end) - 2*1000/f_rede);
idx_reg  = t_ms >= t_regime;

subplot(2,1,1)   % tensões de fase
hold on
plot(t_ms(idx_reg), v_ag(idx_reg), 'b', 'LineWidth', 0.7)
plot(t_ms(idx_reg), v_bg(idx_reg), 'r', 'LineWidth', 0.7)
plot(t_ms(idx_reg), v_cg(idx_reg), 'g', 'LineWidth', 0.7)
% Referências sobrepostas
plot(t_ctrl_ms, v_ref(:,1), 'b--', 'LineWidth', 1.6)
plot(t_ctrl_ms, v_ref(:,2), 'r--', 'LineWidth', 1.6)
plot(t_ctrl_ms, v_ref(:,3), 'g--', 'LineWidth', 1.6)
yline( Vdc/2,  'k:', sprintf('+V_{dc}/2 = +%gV', Vdc/2))
yline(-Vdc/2,  'k:', sprintf('-V_{dc}/2 = -%gV', Vdc/2))
ylabel('Tensão (V)'); grid on
legend('v_{ag}','v_{bg}','v_{cg}','v_a^* (GPC)','v_b^* (GPC)','v_c^* (GPC)', ...
       'Location','best','FontSize',7)
title(sprintf('Tensões de Fase | 2 Níveis: v_{Xg} \\in {±V_{dc}/3, ±2V_{dc}/3}'))
xlim([t_regime t_ms(end)])

subplot(2,1,2)   % tensões de linha
hold on
plot(t_ms(idx_reg), v_ab(idx_reg), 'b', 'LineWidth', 0.7)
plot(t_ms(idx_reg), v_bc(idx_reg), 'r', 'LineWidth', 0.7)
plot(t_ms(idx_reg), v_ca(idx_reg), 'g', 'LineWidth', 0.7)
yline( Vdc, 'k:', sprintf('+V_{dc} = +%gV',Vdc))
yline(-Vdc, 'k:', sprintf('-V_{dc} = -%gV',Vdc))
ylabel('Tensão (V)'); xlabel('Tempo (ms)'); grid on
legend('v_{ab}','v_{bc}','v_{ca}','Location','best','FontSize',8)
title(sprintf('Tensões de Linha | 2 Níveis: v_{LL} \\in {-V_{dc}, 0, +V_{dc}}'))
xlim([t_regime t_ms(end)])

%-- Plot 4: Duty Médio vs Índice de Modulação -----------------------------
figure('Name','[PWM] Duty Médio vs Ma','Color','w', ...
       'Position',[150 300 900 400]);

subplot(2,1,1)
hold on
plot(t_ctrl_ms, duty_media_a, 'b-o', 'MarkerSize',3, 'LineWidth',1.2)
plot(t_ctrl_ms, duty_media_b, 'r-o', 'MarkerSize',3, 'LineWidth',1.2)
plot(t_ctrl_ms, duty_media_c, 'g-o', 'MarkerSize',3, 'LineWidth',1.2)
yline(0.5,'k--','d=0.5 (V_{dc}/2)','LabelHorizontalAlignment','left')
ylim([0 1]); ylabel('Duty médio'); grid on
legend('d_a','d_b','d_c','Location','best')
title('Duty Cycle Médio por Período de Controle')

subplot(2,1,2)
hold on
plot(t_ctrl_ms, abs(Ma_a(1:ratio:end)), 'b-', 'LineWidth',1.2)
plot(t_ctrl_ms, abs(Ma_b(1:ratio:end)), 'r-', 'LineWidth',1.2)
plot(t_ctrl_ms, abs(Ma_c(1:ratio:end)), 'g-', 'LineWidth',1.2)
yline(1.0,'k--','Ma=1 (limite linear)','LabelHorizontalAlignment','left')
ylim([0 1.3]); ylabel('|Ma|'); xlabel('Tempo (ms)'); grid on
legend('|Ma_a|','|Ma_b|','|Ma_c|','Location','best')
title('Índice de Modulação |Ma| = |v_{ref}| / (V_{dc}/2)')

% Relatório no console
fprintf('\n=== applyPWM — 2 Níveis ===\n');
fprintf('  Janela   : k=%d → k=%d (%d instantes de controle)\n', k_inicio, k_fim, n_janela);
fprintf('  t_pwm    : %.3f → %.3f ms (%d amostras)\n', t_ms(1), t_ms(end), n_pwm);
fprintf('  ratio    : %d amostras PWM / ciclo de controle\n', ratio);
fprintf('  Ma_max   : %.4f %s\n', Ma_max, ternary(Ma_max<=1,'[LINEAR]','[SOBREMODULACAO]'));
fprintf('  v_fase   : {0, ±%.1f, ±%.1f} V\n', Vdc/3, 2*Vdc/3);
fprintf('  v_linha  : {0, ±%.1f} V\n', Vdc);
fprintf('===========================\n\n');

end

% Função auxiliar inline
function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
```

---

## Tabela dos 8 Estados de Chaveamento — 2 Níveis
```
Estado | duty_a | duty_b | duty_c | v_ag        | v_bg        | v_cg
-------|--------|--------|--------|-------------|-------------|-------------
  0    |   0    |   0    |   0    |     0       |     0       |     0
  1    |   1    |   0    |   0    |  +2Vdc/3    |  -Vdc/3     |  -Vdc/3
  2    |   0    |   1    |   0    |  -Vdc/3     |  +2Vdc/3    |  -Vdc/3
  3    |   1    |   1    |   0    |  +Vdc/3     |  +Vdc/3     |  -2Vdc/3
  4    |   0    |   0    |   1    |  -Vdc/3     |  -Vdc/3     |  +2Vdc/3
  5    |   1    |   0    |   1    |  +Vdc/3     |  -2Vdc/3    |  +Vdc/3
  6    |   0    |   1    |   1    |  -2Vdc/3    |  +Vdc/3     |  +Vdc/3
  7    |   1    |   1    |   1    |     0       |     0       |     0
