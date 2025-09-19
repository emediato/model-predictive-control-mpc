clear, clc, close all

% Análise de equilíbrio
theta1bar = pi/4;

theta2bar = -pi/6;

[~,~,~,~,tau1bar,tau2bar] = model_fun(theta1bar,theta2bar);


%%% Simula��o com o modelo não linear em malha aberta

% definindo uma chamada para essa function no codigo principal!! 


equacao_estado = @(t,xnl,unl) equacao_estado_function(t, xnl, unl);
% 
xnlbar = [theta1bar;theta2bar;0;0];
unlbar = [tau1bar;tau2bar];

% Linearização e projeto de controlador

%

[Ac,Bc] = model_lin(theta1bar,theta2bar);
T = 0.01; % periodo amostragem
[A,B] = c2dm(Ac,Bc,[],[],T,'zoh');
C = [1 0 0 0;0 1 0 0];
eigc_des = [-20+10j;-20-10j;-5+2j;-5-2j]; % autovalores de MF 
eig_des = exp(eigc_des*T);

K = place(A,B,eig_des);

aux = inv([A - eye(4) B;C zeros(2,2)])*[zeros(4,2);eye(2)];

Nx = aux(1:4,:);
Nu = aux(5:6,:);

rbar = deg2rad([10;10]);
ubar = Nu*rbar
xbar = Nx*rbar


xnl{1} = xnlbar;
trec = 0; xnlrec = xnl{1}';


for k = 1:20
    x{k} = xnl{k} - xnlbar;
    u{k} = - K*(x{k} - xbar) + ubar ;
    unl{k} = u{k} + unlbar;
    
    [tout,  xnlout] = ode45(@(t,xnl) ...
        equacao_estado( t,xnl,unl{k}),[(k-1)*T k*T], xnl{k}); 
        
    %[tout,  xnlout] = dlsim(xnl,unl{k},[(k-1)*T k*T], xnl{k});
     
    trec = [trec;  tout(2:end)];
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
% export_fig -transparent -painters Serie3_Theta1.pdf

figure(2)
h = plot(trec,rad2deg(xnlrec(:,2)));hold on
set(h,'LineWidth',2)
grid
set(gca,'FontSize',13,'GridLineStyle','--')
xlabel('t / s','FontSize',13)
ylabel('\theta_2(t) / graus ','FontSize',13)
% export_fig -transparent -painters Serie3_Theta2.pdf

Unl = [unl{1,:}];
k = 0:200;
figure(3)
h = stairs(k*T,[Unl(1,:) Unl(1,end)]); hold on
set(h,'LineWidth',2)
grid
set(gca,'FontSize',13,'GridLineStyle','--')
xlabel('t / s','FontSize',13)
ylabel('\tau_1(t) / Nm','FontSize',13)
% export_fig -transparent -painters Serie3_Tau1.pdf

figure(4)
h = stairs(k*T,[Unl(2,:) Unl(2,end)]); hold on
set(h,'LineWidth',2)
grid
set(gca,'FontSize',13,'GridLineStyle','--')
xlabel('t / s','FontSize',13)
ylabel('\tau_2(t) / Nm','FontSize',13)
% export_fig -transparent -painters Serie3_Tau2.pdf
