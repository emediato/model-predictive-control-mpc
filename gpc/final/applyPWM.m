%% Função para após calcular sinal de controle, 
%% aplicar PWM por alguns ciclos

function [a] = applyPWM(saidas, nin, nit, Ts, t_control, fsw, Vdc, K_inv, Nu, N1, N2, delta, lambda, K_park)

    % Park (abc → dq) - para medição de correntes
    K_park = @(theta) (2/3) * [cos(theta),   cos(theta - 2*pi/3),   cos(theta - 4*pi/3);
                               -sin(theta), -sin(theta - 2*pi/3), -sin(theta - 4*pi/3)];

    % pré alocação
    t = ((nin:nit)-nin)*Ts;
    ind = nin:nit;
    
    Vdc = 60;            % tensão do barramento
    fsw = 10000;         % freq. PWM (Hz)
    Ts_pwm = 1/fsw;
    Tsim = 0.04;         % tempo de simulação (s)
    
    fs = 200000;         % taxa de amostragem para simulação dos sinais
    t = 0:1/fs:Tsim;
    f0 = 60; omega = 2*pi*f0;

    % Inicializar sinais modulantes
    m_abc = zeros(3, length(t_control)); % [ma; mb; mc]
    ma = zeros(size(t)); mb=ma; mc=ma;
    va_ref = zeros(size(t)); vb_ref=va_ref; vc_ref=va_ref;
    ga = zeros(size(t)); gb=ga; gc=ga;

    fprintf('Vd: %.1f a %.1f V\n', min(entradas(1,ind)), max(entradas(1,ind)));
    fprintf('Vq: %.1f a %.1f V\n', min(entradas(2,ind)), max(entradas(2,ind)));
    
    
    % carrier triangular (centred -1..1)
    %carrier = sawtooth(2*pi*fsw*t, 0.5); % triangular entre -1 e 1
    carrier = 0.5 + 0.5*sawtooth(2*pi*fsw*t, 0.5);  % Triangular [0,1]
    
    u_d_array = saidas(1,ind);
    u_q_array = saidas(2,ind);
    
    % inverse Clarke (αβ -> abc)
    K_inv = [1, 0; -1/2, sqrt(3)/2; -1/2, -sqrt(3)/2];
    
    % pré alocação
    ma = zeros(size(t)); mb=ma; mc=ma;
    va_ref = zeros(size(t)); vb_ref=va_ref; vc_ref=va_ref;
    ga = zeros(size(t)); gb=ga; gc=ga;
    
    idx = 1;
    
    for k=1:ind
        u_d = u_d_array(k);
        u_q = u_q_array(k);
        t_start = (k-1)*Ts_control;
        t_end = min(k*Ts_control, Tsim);
        idxs = find(t>=t_start & t < t_end);
    
        for jj = 1:length(idxs)
            ti = t(idxs(jj));
            theta = omega*ti;                 % ângulo instantâneo (ex.: grid angle)
            % inverse Park
            v_alpha =  cos(theta)*u_d - sin(theta)*u_q;
            v_beta  =  sin(theta)*u_d + cos(theta)*u_q;
    
            % inverse Clarke => va,vb,vc (fase-neutral references)
            vabc = K_inv * [v_alpha; v_beta];
            va = vabc(1); vb = vabc(2); vc = vabc(3);
    
            % modulation indices (SPWM)
            ma_i = 0.5*va / Vdc;
            mb_i = 0.5*vb / Vdc;
            mc_i = 0.5*vc / Vdc;
    
            % limit
            ma_i = max(min(ma_i,1),0);
            mb_i = max(min(mb_i,1),0);
            mc_i = max(min(mc_i,1),0);
    
            % store
            va_ref(idxs(jj)) = va; vb_ref(idxs(jj)) = vb; vc_ref(idxs(jj)) = vc;
            ma(idxs(jj)) = ma_i; mb(idxs(jj)) = mb_i; mc(idxs(jj)) = mc_i;
    
            % compare with carrier to make gate signals
            carrier_val = carrier(idxs(jj));
            ga(idxs(jj)) = ma_i > carrier_val;
            gb(idxs(jj)) = mb_i > carrier_val;
            gc(idxs(jj)) = mc_i > carrier_val;
        end
    end
    
    % plots rápidos
    figure;
    n_plots = 3;
    subplot(n_plots,1,1); 
    plot(t,va_ref,'b', t, vb_ref,'r', t, vc_ref,'g', 'LineWidth',5); 
    ylabel('v^*_ph (V)'); 
    legend('va^*','vb^*','vc^*')
    title(sprintf(['Horizontes: Nu=[%d %d], N1=[%d %d], N2=[%d %d], ', ...
                   'delta=[%.4f %.4f], lambda=[%.4f %.4f]'], ...
                   Nu(1), Nu(2), N1(1), N1(2), N2(1), N2(2), ...
                   delta(1), delta(2), lambda(1), lambda(2)));

    subplot(n_plots,1,2); plot(t,ma,t,mb,t,mc,'LineWidth',5);
    ylabel('mod indices');
    
    subplot(n_plots,1,3); 
    plot(t,ga,'b',t,gb,'r',t,gc,'g','LineWidth',1);
    ylabel('gates (0/1)'); xlabel('time (s)');

    a = 1;
%     subplot(n_plots,1,4); 
%     plot(t,ga,'b',t,gb,'r',t,gc,'g','LineWidth',1);
%     ylabel('gates (0/1)'); xlabel('time (s)')
end
