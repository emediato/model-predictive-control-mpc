
Q = diag([50,1]); R = 1;
Ac = [0 1; 0 0]
Bc=[0; 1]
C=[1 0]
T=0.1
A = [1 T; 0 1]
B=[T^2/2; T]

Sf=0
bf=0

xbar = zeros(4,1); ubar = 0;

[K,P] = dlqr(A,B,Q,R);

un = dlqrcon(A, B, Q, R,P,N, Sx, bx, Sf, bf, x{k}, xbar, ubar);





%% 


N = 15; % Horizonte de controle a ser adotado
x{1} = x0;
sys = ss(Ac,Bc,[],[]);
trec = 0; xrec = x{1}'; 
xbar = zeros(4,1); ubar = 0;
nt = 10; kf = 100; warning off

t = 0:nt:kf; % Vetor de tempo
u = zeros(1, kf); % 1 sinal de controle, cada coluna representa um tempo
y = zeros(1, kf); 


% u(k) = un(1);

% considerando eq de saída da forma y(t)=C.x(t)

C = [0 1 0 0];

k=1;
while (max(Sf*x{1} - bf) > 0)
    un = dlqrcon(A, B, Q, R,P,N, Sx, bx, Sf, bf, x{k}, xbar, ubar);
    u(k) = un(1);
    x{k+1} = A*x{x} + B*u(k);
    k=k+1;
end

kf = k;

for k = 1:kf-1  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Insira linhas de código aqui para obter u{k} com base em x{k} %

    u(k) = ubar - K*(x{k} - xbar) 
    x{k+1} = A*x{k} + B*u(k);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = linspace((k-1)*T,k*T,nt);
    %U = repmat(u{k}',nt,1);
    U = repmat(u(k)',nt,1);
    [~,tout,xout] = lsim(sys,U,t,x{k});
    trec = [trec;tout(2:end)];
    xrec = [xrec;xout(2:end,:)];
    x{k+1} = xrec(end,:)';
    
    y(k) = C * x(:, k);
end
