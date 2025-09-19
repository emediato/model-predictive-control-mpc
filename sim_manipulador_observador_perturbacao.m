clear, clc, close all

% Análise de equilíbrio
theta1bar = pi/4;
theta2bar = -pi/6;
[~,~,~,~,tau1bar,tau2bar] = model_fun(theta1bar,theta2bar);

% %% Simulação com o modelo não linear em malha aberta
equacao_estado = @(t,xnl,unl) equacao_estado_function(t,xnl,unl);
% 
xnlbar = [theta1bar;theta2bar;0;0];
unlbar = [tau1bar;tau2bar];

% Linearização e projeto de controlador

[Ac,Bc] = model_lin(theta1bar,theta2bar);
T = 0.01;
[A,B] = c2dm(Ac,Bc,[],[],T,'zoh');
C = [1 0 0 0;0 1 0 0];
eigc_des = [-20+10j;-20-10j;-5+2j;-5-2j];
eig_des = exp(eigc_des*T);
K = place(A,B,eig_des)

aux = inv([A - eye(4) B;C zeros(2,2)])*[zeros(4,2);eye(2)];
Nx = aux(1:4,:);
Nu = aux(5:6,:);

rbar = deg2rad([10;10]);

% Projeto do observador deadbeat
E = B;
Achi = [A E;zeros(2,4) eye(2)];
Bchi = [B;zeros(2,2)];
H = eye(4);
Hchi = [H zeros(4,2)];
eig_des = [0:5]*1e-2;
L = (place(Achi',Hchi',eig_des))'
M = inv(Achi)*L
eig(Achi - M*Hchi*Achi)

% Relações de equilíbrio considerando a perturbação
aux = inv([A - eye(4) B;C zeros(2,2)])*[-E;zeros(2,2)];
%[6x2]
Mx = aux(1:4,:)
Mu = aux(5:6,:);
E = B;
Achi = [A E;zeros(2,4) eye(2)];
Bchi = [B;zeros(2,2)];
H = eye(4);
Hchi = [H zeros(4,2)];
eig_des = [0:5]*1e-2;
L = (place(Achi',Hchi',eig_des))';
M = inv(Achi)*L;
eig(Achi - M*Hchi*Achi) % 6 AUTOVALORES
aux = inv([A - eye(4) B;C zeros(2,2)])*[-E;zeros(2,2)];
Mx = aux(1:4,:)
Mu = aux(5:6,:);
xnl{1} = xnlbar;
trec = 0; xnlrec = xnl{1}';
chi_pri{1} = zeros(6,1);



% Com observador
clear xbar
clear ubar
xnl{1} = xnlbar;
trec = 0; xnlrec = xnl{1}';
chi_pri{1} = zeros(6,1);
for k = 1:200
    x{k} = xnl{k} - xnlbar;
    zpri{k} = Hchi*chi_pri{k};
    z{k} = x{k};
    chi_pos{k} = chi_pri{k} + M*(z{k} - zpri{k});
    dbar_pos{k} = chi_pos{k}(5:6);
    xbar{k} = Nx*rbar + Mx*dbar_pos{k};
    ubar{k} = Nu*rbar + Mu*dbar_pos{k};
    u{k} = ubar{k} - K*(x{k} - xbar{k});
    chi_pri{k+1} = Achi*chi_pos{k} + Bchi*u{k};
    unl{k} = u{k} + unlbar;
    [tout,xnlout] = ode45(@(t,xnl) equacao_estado(t,xnl,unl{k}),[(k-1)*T k*T],xnl{k}); 
    trec = [trec;tout(2:end)];
    xnlrec = [xnlrec;xnlout(2:end,:)];
    xnl{k+1} = xnlrec(end,:)';
        
end

figure(1)
h = plot(trec,rad2deg(xnlrec(:,1)));hold on
set(h,'LineWidth',2)
grid
set(gca,'FontSize',13,'GridLineStyle','--')
xlabel('t / s','FontSize',13)
ylabel('\theta_1(t) / graus ','FontSize',13)
% export_fig -transparent -painters Serie4_Theta1.pdf

figure(2)
h = plot(trec,rad2deg(xnlrec(:,2)));hold on
set(h,'LineWidth',2)
grid
set(gca,'FontSize',13,'GridLineStyle','--')
xlabel('t / s','FontSize',13)
ylabel('\theta_2(t) / graus ','FontSize',13)
% export_fig -transparent -painters Serie4_Theta2.pdf

Unl = [unl{1,:}];
k = 0:200;
figure(3)
h = stairs(k*T,[Unl(1,:) Unl(1,end)]); hold on
set(h,'LineWidth',2)
grid
set(gca,'FontSize',13,'GridLineStyle','--')
xlabel('t / s','FontSize',13)
ylabel('\tau_1(t) / Nm','FontSize',13)
% export_fig -transparent -painters Serie4_Tau1.pdf

figure(4)
h = stairs(k*T,[Unl(2,:) Unl(2,end)]); hold on
set(h,'LineWidth',2)
grid
set(gca,'FontSize',13,'GridLineStyle','--')
xlabel('t / s','FontSize',13)
ylabel('\tau_2(t) / Nm','FontSize',13)
% export_fig -transparent -painters Serie4_Tau2.pdf
