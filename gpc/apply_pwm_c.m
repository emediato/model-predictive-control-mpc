function applyPWM(in, out, i_init, i_sim, Ts_pwm, Ts_control, ...
                  t_control, fs, Vdc, Nu, N1, N2, delta, lambda)
% =========================================================================
% applyPWM — Modulação PWM Trifásica com entrada do controlador MPC
%
% DESCRIÇÃO DOS PARÂMETROS (conforme dissertação de mestrado):
%
%   in         : struct — sinal de controle ótimo calculado pelo MPC
%                  in.u  [n_ctrl x 3] tensões de referência v_a*, v_b*, v_c*
%                         OU
%                  in.vd, in.vq  referências no referencial dq
%
%   out        : struct — saída da simulação do processo
%                  out.y  [n_ctrl x n_saidas] correntes/tensões medidas
%
%   i_init     : iteração inicial — índice de início (inteiro >= 1)
%                Garante inicialização correta do buffer circular do MPC
%                Importante quando applyPWM é chamada em subjanelas
%
%   i_sim      : número total de iterações da simulação (inteiro)
%                Define o comprimento do loop principal
%
%   Ts_pwm     : período de chaveamento (s) = 1/fs
%   Ts_control : período do controlador  (s) — múltiplo inteiro de Ts_pwm
%
%   t_control  : vetor de instantes de controle [1 x n_ctrl] (s)
%                Gerado externamente: t_control = 0:Ts_control:t_final
%
%   fs         : frequência de chaveamento (Hz)
%   Vdc        : tensão do barramento CC (V)
%
%   Nu         : horizonte de controle do MPC
%                Número de incrementos Δu calculados a cada instante k
%                Nu <= N2-N1+1 (restrição padrão MPC)
%
%   N1         : horizonte de predição inicial (inclui atraso puro d)
%                N1 >= d+1, onde d = floor(atraso_total / Ts_control)
%                Dimensão: escalar (aplicado uniformemente às 3 fases)
%
%   N2         : horizonte de predição final
%                Janela de predição: [N1, N2] — comprimento N2-N1+1
%
%   delta      : ponderação do incremento de controle Δu (escalar >= 0)
%                Alto delta → controle suave, resposta lenta
%                Baixo delta → controle agressivo, possível instabilidade
%
%   lambda     : ponderação dos erros futuros (escalar > 0)
%                Peso na função custo J = lambda*||e||² + delta*||Δu||²
%
% SAÍDAS: plots gerados internamente (dissertação — análise visual)
% =========================================================================

% -------------------------------------------------------------------------
% PASSO 1 — VALIDAÇÃO E CONFIGURAÇÃO INICIAL
% -------------------------------------------------------------------------
% Verificação de consistência entre parâmetros MPC e PWM.
% Erros aqui indicam configuração incorreta do experimento.

% Razão entre períodos: quantos ciclos PWM por ciclo de controle
ratio = round(Ts_control / Ts_pwm);

% Verificação: Ts_control deve ser múltiplo inteiro de Ts_pwm
% (sincronismo entre controlador digital e modulador)
if abs(ratio - Ts_control/Ts_pwm) > 1e-9
    warning('applyPWM: Ts_control não é múltiplo exato de Ts_pwm. Verifique parametrização.');
end

% Verificação do horizonte MPC
if Nu > (N2 - N1 + 1)
    error('applyPWM: Nu > N2-N1+1 viola restrição do MPC. Ajuste horizontes.');
end

% Número de instantes de controle efetivos
n_ctrl = length(t_control);

% Índice de início do loop (i_init permite hot-start da simulação)
% i_init = 1: simulação desde o início
% i_init > 1: continuação de simulação anterior (janela deslizante)
k_start = max(1, i_init);
k_end   = min(i_sim, n_ctrl);   % respeita número de iterações solicitado

fprintf('=== applyPWM iniciada ===\n');
fprintf('  Iterações : %d → %d  (total: %d)\n', k_start, k_end, k_end-k_start+1);
fprintf('  fs         = %.1f kHz\n',  fs/1e3);
fprintf('  Ts_pwm     = %.2f µs\n',   Ts_pwm*1e6);
fprintf('  Ts_control = %.2f µs\n',   Ts_control*1e6);
fprintf('  ratio      = %d PWM/ctrl\n', ratio);
fprintf('  Horizonte  : N1=%d, N2=%d, Nu=%d\n', N1, N2, Nu);
fprintf('  Pesos MPC  : lambda=%.4f, delta=%.4f\n', lambda, delta);

% -------------------------------------------------------------------------
% PASSO 2 — EXTRAÇÃO DOS SINAIS DE CONTROLE ÓTIMO (in)
% -------------------------------------------------------------------------
% "in" contém o sinal de controle ótimo calculado pelo MPC.
% O MPC resolve min J = lambda*sum(e²) + delta*sum(Δu²)
% e retorna a sequência de tensões de referência v_abc*(k).
%
% Convenção adotada: in.u é matriz [n_ctrl x 3] com colunas [va* vb* vc*]
% Se disponível in.vd e in.vq, converte para abc via Park⁻¹ + Clarke⁻¹

if isfield(in, 'u') && size(in.u, 2) >= 3
    % Caso 1: MPC fornece diretamente v_abc* (referências de fase)
    v_ref_abc = in.u(k_start:k_end, 1:3);   % [n_iter x 3]

elseif isfield(in, 'vd') && isfield(in, 'vq')
    % Caso 2: MPC fornece v_dq* → converte para abc
    % Útil quando o MPC opera no referencial síncrono dq
    we_grid = 2*pi*50;   % frequência angular da rede (rad/s)
    n_iter  = k_end - k_start + 1;
    v_ref_abc = zeros(n_iter, 3);

    for k = 1:n_iter
        k_abs   = k + k_start - 1;
        theta_e = we_grid * t_control(k_abs);   % ângulo elétrico

        % Transformada de Park Inversa: dq → αβ
        % [vα]   [ cos θ  -sin θ] [vd]
        % [vβ] = [ sin θ   cos θ] [vq]
        vd = in.vd(k_abs);
        vq = in.vq(k_abs);
        v_alpha =  vd*cos(theta_e) - vq*sin(theta_e);
        v_beta  =  vd*sin(theta_e) + vq*cos(theta_e);

        % Transformada de Clarke Inversa: αβ → abc (amplitude invariante)
        % va =  vα
        % vb = -vα/2 + (√3/2)·vβ
        % vc = -vα/2 - (√3/2)·vβ
        v_ref_abc(k,1) =  v_alpha;
        v_ref_abc(k,2) = -v_alpha/2 + (sqrt(3)/2)*v_beta;
        v_ref_abc(k,3) = -v_alpha/2 - (sqrt(3)/2)*v_beta;
    end

else
    error('applyPWM: struct "in" deve conter campo "u" [n x 3] ou campos "vd" e "vq".');
end

% -------------------------------------------------------------------------
% PASSO 3 — EXTRAÇÃO DA SAÍDA DO PROCESSO (out)
% -------------------------------------------------------------------------
% "out" contém a saída medida do processo (correntes de fase).
% Usado para: (a) verificação de rastreamento, (b) cálculo de THD de corrente.
%
% Convenção: out.y é matriz [n_ctrl x n_saidas]
% Correntes de fase: colunas 1, 2, 3 → ia, ib, ic

if isfield(out, 'y') && ~isempty(out.y)
    i_medido = out.y(k_start:k_end, :);   % correntes medidas
    tem_corrente = true;
else
    tem_corrente = false;
    warning('applyPWM: out.y não encontrado — plots de corrente omitidos.');
end

% -------------------------------------------------------------------------
% PASSO 4 — CONSTRUÇÃO DO EIXO DE TEMPO PWM
% -------------------------------------------------------------------------
% O eixo PWM tem resolução Ts_pwm (muito mais fino que Ts_control).
% Cada intervalo de controle [k, k+1] contém "ratio" períodos PWM.
%
% t_pwm é construído para as iterações k_start → k_end.

t_inicio = t_control(k_start);
t_fim    = t_control(k_end);
t_pwm    = t_inicio : Ts_pwm : t_fim;
n_pwm    = length(t_pwm);

% -------------------------------------------------------------------------
% PASSO 5 — PORTADORA TRIANGULAR (SPWM Natural Sampling)
% -------------------------------------------------------------------------
% A portadora triangular simétrica é o padrão industrial para SPWM.
% Características:
%   - Frequência = fs (chaveamento)
%   - Amplitude  = ±1 (normalizada — referência será normalizada por Vdc/2)
%   - Fase       = 0 (referência arbitrária)
%
% sawtooth(t, 0.5): dente-de-serra com pico em 50% do período → triângulo
%
% Justificativa triangular simétrica:
%   - Cancela harmônicos de ordem par (simetria de meia-onda)
%   - Concentra harmônicos em torno de fs e múltiplos (2fs, 3fs...)
%   - Harmônicos distantes da fundamental → fácil filtragem

carrier = sawtooth(2*pi*fs*t_pwm, 0.5);
% Resultado: vetor [1 x n_pwm] com valores em [-1, +1]

% -------------------------------------------------------------------------
% PASSO 6 — INTERPOLAÇÃO ZOH (Zero-Order Hold)
% -------------------------------------------------------------------------
% As referências v_abc* existem apenas nos instantes de controle (Ts_control).
% O modulador PWM opera em Ts_pwm (mais rápido).
%
% ZOH ('previous'): mantém o valor constante até a próxima atualização.
% Justificativa: emula o comportamento real do DSP — o registrador de
% comparação do timer PWM é atualizado apenas quando o controlador executa.
% Qualquer outra interpolação (linear, spline) seria fisicamente incorreta.

t_ctrl_local = t_control(k_start:k_end);
n_iter_local = length(t_ctrl_local);

% Expansão ZOH: cada valor de controle mantido por "ratio" amostras PWM
va_pwm = interp1(t_ctrl_local, v_ref_abc(:,1), t_pwm, 'previous', 'extrap');
vb_pwm = interp1(t_ctrl_local, v_ref_abc(:,2), t_pwm, 'previous', 'extrap');
vc_pwm = interp1(t_ctrl_local, v_ref_abc(:,3), t_pwm, 'previous', 'extrap');

% -------------------------------------------------------------------------
% PASSO 7 — NORMALIZAÇÃO E ÍNDICE DE MODULAÇÃO
% -------------------------------------------------------------------------
% O comparador PWM compara a referência normalizada com a portadora [-1,+1].
% Normalização: Ma = v_ref / (Vdc/2)
%
% Limite de modulação linear: Ma ≤ 1
%   Para Ma > 1 → sobremodulação → harmônicos de baixa ordem indesejados
%
% Verificação do índice de modulação (IM):
%   Para SPWM trifásico: tensão de fase máxima = Vdc/2
%   Para SVPWM:          tensão de fase máxima = Vdc/√3 ≈ 0.577*Vdc

Ma_a = va_pwm / (Vdc/2);
Ma_b = vb_pwm / (Vdc/2);
Ma_c = vc_pwm / (Vdc/2);

% Índice de modulação médio (para log)
Ma_rms = sqrt(mean(Ma_a.^2 + Ma_b.^2 + Ma_c.^2)/3);
if max(abs([Ma_a, Ma_b, Ma_c])) > 1
    warning('applyPWM: Sobremodulação detectada (Ma > 1). Verifique Vdc ou referências.');
end

% -------------------------------------------------------------------------
% PASSO 8 — COMPARAÇÃO COM PORTADORA → DUTY CYCLES
% -------------------------------------------------------------------------
% Lógica de chaveamento SPWM:
%   Se Ma(k) >= carrier(k) → chave superior LIGADA → saída = +Vdc/2
%   Se Ma(k) <  carrier(k) → chave inferior LIGADA → saída = -Vdc/2
%
% duty = 1: chave superior conduzindo (polo em +Vdc/2)
% duty = 0: chave inferior conduzindo (polo em -Vdc/2)

duty_a = double(Ma_a >= carrier);   % [1 x n_pwm] binário
duty_b = double(Ma_b >= carrier);
duty_c = double(Ma_c >= carrier);

% -------------------------------------------------------------------------
% PASSO 9 — RECONSTRUÇÃO DAS TENSÕES DE FASE
% -------------------------------------------------------------------------
% Tensão de polo (phase-to-negative DC bus):
%   v_XN = duty_X * Vdc
%   Para duty=1: v_XN = Vdc  |  Para duty=0: v_XN = 0
%
% Tensão de fase (phase-to-neutral, assumindo neutro virtual equilibrado):
%   v_ag = (2*v_AN - v_BN - v_CN) / 3
%
% Justificativa da fórmula: deriva da condição de soma nula (v_ag+v_bg+v_cg=0)
% e das relações v_XN = v_Xg + v_gN para neutro flutuante.

v_AN = duty_a * Vdc;
v_BN = duty_b * Vdc;
v_CN = duty_c * Vdc;

% Tensões de fase (line-to-neutral)
v_ag = (2*v_AN - v_BN - v_CN) / 3;
v_bg = (2*v_BN - v_AN - v_CN) / 3;
v_cg = (2*v_CN - v_AN - v_BN) / 3;

% Tensões de linha (line-to-line)
v_ab = v_AN - v_BN;   % = v_ag - v_bg
v_bc = v_BN - v_CN;
v_ca = v_CN - v_AN;

% -------------------------------------------------------------------------
% PASSO 10 — CÁLCULO DO THD (Total Harmonic Distortion)
% -------------------------------------------------------------------------
% THD é métrica de qualidade obrigatória em dissertações de EP.
% Norma IEC 61000-3-2 limita THD de corrente para cargas < 16A.
% IEEE 519-2022 define limites por nível de tensão.
%
% Método: FFT → componente fundamental → harmônicos até 50ª ordem
%
% Resolução espectral: Δf = fs_amostragem / N_fft
% Para boa resolução em 50 Hz: N_fft deve ser múltiplo de fs/50

N_fft = n_pwm;
df    = 1 / (N_fft * Ts_pwm);   % resolução em frequência (Hz)
f_ax  = (0:N_fft-1) * df;       % eixo de frequências

% FFT das tensões de fase A
Y_ag    = fft(v_ag) / N_fft;
espectro = 2 * abs(Y_ag(1:floor(N_fft/2)+1));
f_meio  = f_ax(1:floor(N_fft/2)+1);

% Identificação da fundamental (50 Hz)
[~, idx_f1] = min(abs(f_meio - 50));
V1 = espectro(idx_f1);           % amplitude da fundamental

% Cálculo THD (até 50ª harmônica — padrão IEC/IEEE)
idx_h_max = min(50*idx_f1, length(espectro));
espectro_harm = espectro;
espectro_harm(idx_f1) = 0;       % zera a fundamental

THD_v = (sqrt(sum(espectro_harm(1:idx_h_max).^2)) / V1) * 100;   % em %

% THD de corrente (se disponível)
THD_i = NaN;
if tem_corrente && size(i_medido,2) >= 1
    ia_ctrl = i_medido(:,1);
    % Interpola corrente para resolução PWM
    ia_pwm = interp1(t_ctrl_local, ia_ctrl, t_pwm, 'linear', 'extrap');

    Yi    = fft(ia_pwm) / n_pwm;
    esp_i = 2 * abs(Yi(1:floor(n_pwm/2)+1));
    I1    = esp_i(idx_f1);
    if I1 > 1e-6
        esp_i_harm = esp_i;
        esp_i_harm(idx_f1) = 0;
        THD_i = (sqrt(sum(esp_i_harm(1:idx_h_max).^2)) / I1) * 100;
    end
end

fprintf('\n  THD tensão fase A = %.2f%%\n', THD_v);
if ~isnan(THD_i)
    fprintf('  THD corrente fase A = %.2f%%\n', THD_i);
end

% =========================================================================
% GERAÇÃO DE PLOTS (dentro da função, conforme solicitado)
% =========================================================================

t_ms  = t_pwm * 1e3;       % eixo de tempo em ms (melhor leitura)
t_ctrl_ms = t_ctrl_local * 1e3;

% Janela de visualização detalhada: últimos 3 ciclos de rede (60 ms)
% Justificativa: primeiros ciclos contêm transitório de partida (MPC)
% A análise em regime permanente é mais relevante para o THD
f_rede    = 50;
t_regime  = max(t_ms(1), t_ms(end) - 3*(1000/f_rede));
idx_zoom  = t_ms >= t_regime;
idx_zoom_ctrl = t_ctrl_ms >= t_regime;

% -----------------------------------------------------------------
% PLOT 1 — Sinais de Referência do MPC (v_abc*)
% -----------------------------------------------------------------
% Mostra o que o otimizador MPC calculou como tensão ideal.
% Deve ser aproximadamente senoidal em regime permanente.
% Distorções indicam saturação ou atraso de predição.

figure('Name','[applyPWM] Referências MPC — v_abc*', ...
       'NumberTitle','off','Color','w','Position',[100 600 900 400]);

subplot(3,1,1)
plot(t_ctrl_ms, v_ref_abc(:,1), 'b-', 'LineWidth', 1.2)
ylabel('v_a^* (V)'); grid on; xlim([t_ctrl_ms(1) t_ctrl_ms(end)])
title('Referências de Tensão calculadas pelo MPC — v_{abc}^*')
xline(t_ctrl_ms(find(idx_zoom_ctrl,1)), 'r--', 'Regime', 'LabelOrientation','horizontal')

subplot(3,1,2)
plot(t_ctrl_ms, v_ref_abc(:,2), 'r-', 'LineWidth', 1.2)
ylabel('v_b^* (V)'); grid on; xlim([t_ctrl_ms(1) t_ctrl_ms(end)])

subplot(3,1,3)
plot(t_ctrl_ms, v_ref_abc(:,3), 'g-', 'LineWidth', 1.2)
ylabel('v_c^* (V)'); xlabel('Tempo (ms)'); grid on
xlim([t_ctrl_ms(1) t_ctrl_ms(end)])

% -----------------------------------------------------------------
% PLOT 2 — Portadora Triangular vs Referências Normalizadas
% -----------------------------------------------------------------
% Visualização central do SPWM: a interseção entre portadora e
% referência determina os instantes de chaveamento.
% Apresentado em zoom (3 períodos PWM) para legibilidade.

n_zoom_pwm  = min(round(5/fs/Ts_pwm), n_pwm);   % 5 períodos PWM
idx_carr    = (n_pwm - n_zoom_pwm + 1) : n_pwm;  % últimos 5 períodos

figure('Name','[applyPWM] Portadora vs Referências (SPWM)', ...
       'NumberTitle','off','Color','w','Position',[100 150 900 500]);

subplot(2,1,1)
t_carr_us = t_pwm(idx_carr) * 1e6;   % µs para zoom
plot(t_carr_us, carrier(idx_carr), 'k-', 'LineWidth', 1.5); hold on
plot(t_carr_us, Ma_a(idx_carr),    'b-', 'LineWidth', 1.8)
plot(t_carr_us, Ma_b(idx_carr),    'r-', 'LineWidth', 1.8)
plot(t_carr_us, Ma_c(idx_carr),    'g-', 'LineWidth', 1.8)
yline( 1, 'k:', 'LineWidth', 0.8)
yline(-1, 'k:', 'LineWidth', 0.8)
ylim([-1.3 1.3])
legend('Portadora','M_a = v_a^*/(V_{dc}/2)', ...
       'M_b','M_c','Location','best')
xlabel('Tempo (µs)'); ylabel('Amplitude Normalizada')
title(sprintf('Portadora Triangular vs Referências — f_s = %g kHz | Zoom: %d períodos PWM', ...
              fs/1e3, 5))
grid on

subplot(2,1,2)
stairs(t_carr_us, duty_a(idx_carr), 'b-', 'LineWidth', 1.5); hold on
stairs(t_carr_us, duty_b(idx_carr)+2, 'r-', 'LineWidth', 1.5)
stairs(t_carr_us, duty_c(idx_carr)+4, 'g-', 'LineWidth', 1.5)
yticks([0 1 2 3 4 5])
yticklabels({'0','1','0','1','0','1'})
ylabel('Duty (a / b / c)'); xlabel('Tempo (µs)')
title('Sinais de Chaveamento — Duty Cycles (empilhados)')
grid on

% -----------------------------------------------------------------
% PLOT 3 — Tensões de Fase Chaveadas (PWM real)
% -----------------------------------------------------------------
% Forma de onda real na saída do inversor — alta frequência de comutação
% visível. Em dissertação: comparar com referência senoidal.

figure('Name','[applyPWM] Tensões de Fase — Saída PWM', ...
       'NumberTitle','off','Color','w','Position',[100 100 1000 600]);

subplot(3,1,1)
plot(t_ms(idx_zoom), v_ag(idx_zoom), 'b', 'LineWidth', 0.6); hold on
plot(t_ctrl_ms(idx_zoom_ctrl), v_ref_abc(idx_zoom_ctrl,1), ...
     'r--', 'LineWidth', 1.5)
ylabel('v_{ag} (V)'); grid on
title(sprintf('Tensões de Fase — Saída do Inversor | THD(v_{ag}) = %.2f%%', THD_v))
legend('v_{ag} PWM','v_a^* (ref MPC)','Location','best')
xlim([t_ms(find(idx_zoom,1)) t_ms(end)])

subplot(3,1,2)
plot(t_ms(idx_zoom), v_bg(idx_zoom), 'r', 'LineWidth', 0.6); hold on
plot(t_ctrl_ms(idx_zoom_ctrl), v_ref_abc(idx_zoom_ctrl,2), ...
     'b--', 'LineWidth', 1.5)
ylabel('v_{bg} (V)'); grid on
legend('v_{bg} PWM','v_b^* (ref MPC)','Location','best')
xlim([t_ms(find(idx_zoom,1)) t_ms(end)])

subplot(3,1,3)
plot(t_ms(idx_zoom), v_cg(idx_zoom), 'g', 'LineWidth', 0.6); hold on
plot(t_ctrl_ms(idx_zoom_ctrl), v_ref_abc(idx_zoom_ctrl,3), ...
     'k--', 'LineWidth', 1.5)
ylabel('v_{cg} (V)'); xlabel('Tempo (ms)'); grid on
legend('v_{cg} PWM','v_c^* (ref MPC)','Location','best')
xlim([t_ms(find(idx_zoom,1)) t_ms(end)])

% -----------------------------------------------------------------
% PLOT 4 — Tensões de Linha (line-to-line)
% -----------------------------------------------------------------
% As tensões de linha mostram o nível multinível (para Nu > 2)
% e são usadas para verificar equilíbrio trifásico.
% Vab_pico_teórico = Vdc * sqrt(3)/sqrt(2) ≈ 0.612*Vdc (SPWM)

figure('Name','[applyPWM] Tensões de Linha', ...
       'NumberTitle','off','Color','w','Position',[150 100 1000 500]);

subplot(3,1,1)
plot(t_ms(idx_zoom), v_ab(idx_zoom), 'b', 'LineWidth', 0.7)
ylabel('v_{ab} (V)'); grid on
title(sprintf('Tensões de Linha — V_{dc} = %g V | f_s = %g kHz', Vdc, fs/1e3))
yline( Vdc, 'r:', 'V_{dc}',  'LabelHorizontalAlignment','left')
yline(-Vdc, 'r:', '-V_{dc}', 'LabelHorizontalAlignment','left')
xlim([t_ms(find(idx_zoom,1)) t_ms(end)])

subplot(3,1,2)
plot(t_ms(idx_zoom), v_bc(idx_zoom), 'r', 'LineWidth', 0.7)
ylabel('v_{bc} (V)'); grid on
xlim([t_ms(find(idx_zoom,1)) t_ms(end)])

subplot(3,1,3)
plot(t_ms(idx_zoom), v_ca(idx_zoom), 'g', 'LineWidth', 0.7)
ylabel('v_{ca} (V)'); xlabel('Tempo (ms)'); grid on
xlim([t_ms(find(idx_zoom,1)) t_ms(end)])

% -----------------------------------------------------------------
% PLOT 5 — Correntes de Saída (se disponíveis em out.y)
% -----------------------------------------------------------------

if tem_corrente && size(i_medido, 2) >= 3
    figure('Name','[applyPWM] Correntes de Saída', ...
           'NumberTitle','off','Color','w','Position',[200 100 1000 500]);

    labels_i = {'i_a','i_b','i_c'};
    cores_i  = {'b','r','g'};

    for fase = 1:3
        subplot(3,1,fase)
        plot(t_ctrl_ms(idx_zoom_ctrl), i_medido(idx_zoom_ctrl,fase), ...
             cores_i{fase}, 'LineWidth', 1.4)
        ylabel(sprintf('%s (A)', labels_i{fase})); grid on
        if fase == 1
            if ~isnan(THD_i)
                title(sprintf('Correntes de Saída (out.y) | THD(i_a) = %.2f%%', THD_i))
            else
                title('Correntes de Saída (out.y)')
            end
        end
        if fase == 3, xlabel('Tempo (ms)'); end
        xlim([t_ctrl_ms(find(idx_zoom_ctrl,1)) t_ctrl_ms(end)])
    end
end

% -----------------------------------------------------------------
% PLOT 6 — Espectro de Frequência (FFT) — v_ag e i_a
% -----------------------------------------------------------------
% Gráfico obrigatório em dissertações de EP para análise harmônica.
% Referência: IEEE 519-2022 e IEC 61000-3-2.
%
% Eixo X em kHz para visualizar harmônicos de chaveamento (fs, 2fs...).
% Escala log no eixo Y revela harmônicos de baixa amplitude.

figure('Name','[applyPWM] Análise Espectral — FFT', ...
       'NumberTitle','off','Color','w','Position',[250 100 1000 600]);

if tem_corrente && ~isnan(THD_i)
    n_sub = 2;
else
    n_sub = 1;
end

subplot(n_sub, 1, 1)
stem(f_meio/1e3, espectro, 'filled', ...
     'MarkerSize', 2, 'Color', [0.1 0.3 0.8], 'LineWidth', 0.8)
xlim([0, min(3*fs/1e3, f_meio(end)/1e3)])
xlabel('Frequência (kHz)'); ylabel('|V_{ag}| (V)')
title(sprintf('Espectro de Tensão — v_{ag} | THD = %.2f%%', THD_v))
xline(50e-3,  'r-', '50 Hz',        'LabelVerticalAlignment','bottom')
xline(fs/1e3, 'g-', sprintf('%gkHz',fs/1e3), 'LabelVerticalAlignment','bottom')
xline(2*fs/1e3,'m-',sprintf('%gkHz',2*fs/1e3),'LabelVerticalAlignment','bottom')
set(gca,'YScale','log'); grid on

if tem_corrente && ~isnan(THD_i)
    Yi_plot   = fft(ia_pwm) / n_pwm;
    esp_i_plt = 2 * abs(Yi_plot(1:floor(n_pwm/2)+1));

    subplot(2,1,2)
    stem(f_meio/1e3, esp_i_plt, 'filled', ...
         'MarkerSize', 2, 'Color', [0.8 0.1 0.1], 'LineWidth', 0.8)
    xlim([0, min(3*fs/1e3, f_meio(end)/1e3)])
    xlabel('Frequência (kHz)'); ylabel('|I_a| (A)')
    title(sprintf('Espectro de Corrente — i_a | THD_i = %.2f%%', THD_i))
    xline(50e-3,  'r-', '50 Hz',        'LabelVerticalAlignment','bottom')
    xline(fs/1e3, 'g-', sprintf('%gkHz',fs/1e3), 'LabelVerticalAlignment','bottom')
    set(gca,'YScale','log'); grid on
end

% -----------------------------------------------------------------
% PLOT 7 — Diagrama de Parâmetros MPC (resumo visual)
% -----------------------------------------------------------------
% Plot informativo para dissertação: mostra a janela de predição MPC
% em relação ao horizonte temporal. Útil para seção de metodologia.

figure('Name','[applyPWM] Parâmetros MPC — Horizontes', ...
       'NumberTitle','off','Color','w','Position',[300 400 700 300]);

hold on; axis off
xlim([0 10]); ylim([0 6])

% Linha do tempo
annotation('arrow',[0.07 0.93],[0.15 0.15],'Color','k','LineWidth',1.5)
text(9.5, 0.3, 'k', 'FontSize',11, 'FontWeight','bold')

% Marcadores dos horizontes
for h = 0:8
    plot([h h]+1, [0.4 0.6],'k-','LineWidth',1.2)
    text(h+1, 0.1, num2str(h), 'HorizontalAlignment','center','FontSize',8)
end

% Bloco Nu (horizonte de controle)
rectangle('Position',[1, 1.5, Nu, 1], ...
          'FaceColor',[0.2 0.6 1 0.5],'EdgeColor','b','LineWidth',1.5)
text(1 + Nu/2, 2.05, sprintf('N_u = %d\n(Controle)', Nu), ...
     'HorizontalAlignment','center','FontWeight','bold','Color','b','FontSize',9)

% Bloco N1→N2 (horizonte de predição)
rectangle('Position',[N1, 3.0, N2-N1+1, 1], ...
          'FaceColor',[1 0.6 0.2 0.5],'EdgeColor',[0.8 0.4 0],'LineWidth',1.5)
text(N1 + (N2-N1+1)/2, 3.55, ...
     sprintf('N_1=%d → N_2=%d\n(Predição)', N1, N2), ...
     'HorizontalAlignment','center','FontWeight','bold', ...
     'Color',[0.7 0.3 0],'FontSize',9)

% Anotações de pesos
text(0.5, 4.8, sprintf('\\lambda = %.3f  (peso erro futuro)', lambda), ...
     'FontSize', 9, 'Color', [0.7 0.3 0])
text(0.5, 4.2, sprintf('\\delta = %.3f  (peso \\Deltau)', delta), ...
     'FontSize', 9, 'Color', [0.1 0.4 0.9])

title(sprintf('Configuração MPC — i_{init}=%d | i_{sim}=%d iterações', ...
              i_init, i_sim), 'FontSize', 10)
 
fprintf('\n=== Relatório Final applyPWM ===\n');
fprintf('  Iterações executadas : %d\n',    k_end - k_start + 1);
fprintf('  Pontos PWM gerados   : %d\n',    n_pwm);
fprintf('  Índice Ma (RMS)      : %.4f\n',  Ma_rms);
fprintf('  THD tensão (v_ag)    : %.2f%%\n', THD_v);
if ~isnan(THD_i)
fprintf('  THD corrente (i_a)   : %.2f%%\n', THD_i);
end
fprintf('  Plots gerados        : 7\n');
fprintf('================================\n\n');

end
