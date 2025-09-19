Q = diag([50,1]); R = 1;
Ac = [0 1; 0 0]
Bc=[0; 1]
C=[1 0]
T=0.1
A = [1 T; 0 1]
B=[T^2/2; T]

[K,P] = dlqr(A,B,Q,R);

q=1
% P = SOL DA EQ RICATI
rbar = inputs(1:q); % Considera-se rbar = r(kT)
N = 12; % horizonte de tempo/predição/controle

%%

Af = A - B*K

n = 2; q = 1; p = 1;
Achif = [Af B*(K*Nx + Nu);
zeros(q,n) eye(q)];
Gamma = [eye(n) zeros(n,q);
-K (K*Nx + Nu);
zeros(n,n) Nx;
zeros(p,n) Nu];
Spsi = blkdiag(Sx,Su,Sx,Su);
epsilonx = 1e-3*ones(size(Sx,1),1);
epsilonu = 1e-3*ones(size(Su,1),1);
bpsi = [bx;bu;bx-epsilonx;bu-epsilonu];
n = size(B,1); nd = size(E,2); q = size(So,2) - n - nd;

%%

% passo - a cada par de 10 intervalos temos o passo 0.01
x{1} = [0;0];
dt = 0.01;
nt = (T/dt) + 1;
trec = 0; xrec = x{1}';
for k = 1:40
    un = otimizador_u(A,B,Q,R,Pf,N, Sx,bx,Su,bu,Sf,bf,x{k},xbar,ubar);
    u(k) = un(1);
    t = linspace((k-1)*T,k*T,nt)';
    U = repmat(u(k),nt,1);
    [~,~,xout] = lsim(sys,U,t,x{k});
    trec = [trec;t(2:end)];
    xrec = [xrec;xout(2:end,:)];
    x{k+1} = xrec(end,:)';
end
