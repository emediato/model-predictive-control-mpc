diagrama fluxo


A. A. Ahmed, B. K. Koh, and Y. I. Lee, “A comparison of finite control
set and continuous control set model predictive control schemes for speed
control of induction motors,” IEEE Trans. Ind. Informat., vol. 14, no. 4,
pp. 1334–1346, Apr. 2018, doi: 10.1109/TII.2017.2758393.


## IMPLEMENTAÇÃO PRÁTICA NO MATLAB

```matlab
%% IMPLEMENTAÇÃO PRÁTICA NO MATLAB
fprintf('\n IMPLEMENTAÇÃO PRÁTICA \n\n');

% Configuração de tempo
f_pwm = 10000;
f_control = 1000;
Ts_pwm = 1/f_pwm;
ratio = f_pwm / f_control;

% Tempos de processamento (estimativas)
t_calc = 0.0002;    % 200 μs
t_measure = 0.0001; % 100 μs

% Simulação de 2 ciclos completos
fprintf('SIMULAÇÃO DE 2 CICLOS COMPLETOS:\n\n');

ciclo = 1;
tempo_atual = 0;

while ciclo <= 2
    fprintf('--- CICLO %d ---\n', ciclo);
    
    % Fase 1: Cálculo do GPC
    fprintf('t = %.3f ms: Calculando GPC...\n', tempo_atual*1000);
    tempo_atual = tempo_atual + t_calc;
    fprintf('t = %.3f ms: GPC calculado → u = [Vd; Vq]\n', tempo_atual*1000);
    
    % Fase 2: Aplicação do PWM por vários ciclos
    fprintf('t = %.3f ms: Iniciando PWM (%d ciclos)...\n', tempo_atual*1000, ratio);
    
    for pwm_cycle = 1:ratio
        tempo_atual = tempo_atual + Ts_pwm;
        fprintf('t = %.3f ms: PWM ciclo %d/%d\n', tempo_atual*1000, pwm_cycle, ratio);
        
        % Aqui seria gerado o PWM real
        % generate_pwm(u_current);
    end
    
    fprintf('t = %.3f ms: PWM finalizado\n', tempo_atual*1000);
    
    % Fase 3: Medição
    fprintf('t = %.3f ms: Medindo correntes...\n', tempo_atual*1000);
    tempo_atual = tempo_atual + t_measure;
    fprintf('t = %.3f ms: Medição concluída\n', tempo_atual*1000);
    
    ciclo = ciclo + 1;
    fprintf('Tempo parcial: %.3f ms\n\n', tempo_atual*1000);
end

fprintf('=== TEMPO TOTAL DECORRIDO: %.3f ms ===\n', tempo_atual*1000);
