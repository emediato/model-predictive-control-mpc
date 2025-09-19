m1 = 5; m2 = 3; ks = 10; g = 9.8;
Ac = [0 0 1 0;
      0 0 0 1;
      -ks/m1 ks/m1 0 0;
      ks/m2 -ks/m2 0 0];
Bc = [0;0;1/m1;0];
Ec = [0;0;-1;-1];
C = [0 1 0 0];
T = 0.2; % sample


dbar = g* sin(x2) % perturbação cte


% E=Ec * dbar %%%%%%%%%%%%%%%%%%%%%%%%% verify
Sysc = ss(Ac, Bc, [], []);
Sys=c2d(Sysc, T, 'zoh');
pppp
A=Sys.A
B=Sys.B

[~, E] =c2dm(Ac, Ec, [], [], T, 'zoh')

aux = inv([A - eye(4) B; C 0])*[0; 0; 0; 0; 1]

Nx = aux(1:4)
Nu = aux(5)

rbar=1e-1;
ubar=Nu*rbar % + M*dbar ;
xbar=Nx*rbar % + Mx*dbar ;

q=1 % Q = Q_BAR and R = r * R_BAR


x1 = 1
x2 = 1
x3_max = 0.8/6 % [m/s]
x4_max = x3_max
x1_x2_max = 0.05
u_max = 2 % rho


u1 = 1/x1^2
u2 = u1
u3 = 1/(x3_max^2)
u4 = 1/(x4_max^2)
u5 = 1/(x1_x2_max^2) % spring noise affects u1 and u2
rho= 1/(u_max^2)

R= rho
Q = [u1+u5 -u5 0 0; -u5 u2+u5 0 0;0 0 u3 0; 0 0 0 u4]


[K, S, e] = lqr(A,B,Q,R)
[K, S, e] = dlqr(A,B,Q,R)

g = 9.8;
gamma = 0.5;
dbar = g*sind(gamma);

nRow=4;
u = cell(nRow,1) ;
x = cell(nRow,1) ;

u{1}=[0;0;0;0]; 
x{1}=[0;0;0;0]; 

t=[0:Ts:100];


sys = ss(Ac,Bc,C,0);
dt = T/10;
nt = (T/dt) + 1;
x{1} = [0;0];
trec = 0; xrec = x{1}'; yrec = 0;
for k = 1:20
    u{k} = ubar - K*(x{k} - xbar);
    t = linspace((k-1)*T,k*T,nt);
    U = repmat(u{k}',nt,1);
    [yout,tout,xout] = lsim(sys,U,t,x{k});
    trec = [trec;tout(2:end)];
    xrec = [xrec;xout(2:end,:)];
    yrec = [yrec;yout(2:end)];
    x{k+1} = xrec(end,:)';
end



figure(1)


plot(k*T, X(1,:));
% plot(X(1,:),X(2,:),'b-o'); hold on
%x = x(1,:);
%h = plot(k*T, x) ; hold on
%set(h,'LineWidth',2)
grid
set(gca,'FontSize',13,'GridLineStyle','--')
xlabel('t / s','FontSize',13)
ylabel('\x(t) / graus ','FontSize',13)
%export_fig -transparent -painters Serie1_Theta1.pdf
%
figure(2)
h = plot(k*T, X(2,:)); hold on
set(h,'LineWidth',2)
grid
set(gca,'FontSize',13,'GridLineStyle','--')
xlabel('t / s','FontSize',13)
ylabel('\theta_2(t) / graus ','FontSize',13)
% export_fig -transparent -painters Serie3_Theta2.pdf


Unl = [u{1,:}];
k = 0:100;
figure(3)
h = plot(k*T, Unl(1,:)); hold on
set(h,'LineWidth',2)
grid
set(gca,'FontSize',13,'GridLineStyle','--')
xlabel('t / s','FontSize',13)
ylabel('\u(t) ','FontSize',13)
% export_fig -transparent -painters Serie3_Tau1.pdf


