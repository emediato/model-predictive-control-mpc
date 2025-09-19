
Q = diag([50,1]); R = 1;
Ac = [0 1; 0 0]
Bc=[0; 1]
C=[1 0]
T=0.1
A = [1 T; 0 1]
B=[T^2/2; T]

Sx = [1 -1 0 0;-1 1 0 0];
bx = [0.05;0.05];
Sf=[1 1]
bf=[1 0]




[K,P] = dlqr(A,B,Q,R);
x{1} = [0; 0];

N = 15; % Horizonte de controle a ser adotado
%xbar = zeros(4,1); ubar = 0;
xbar = zeros(4,1); ubar = 0;

rbar=1;
un = dlqrcon(A, B, Q, R,P,N, Sx, bx, Sf, bf, x{k}, xbar, ubar);


while (max(Sf*x{1} - bf) > 0)
    un = dlqrcon(A, B, Q, R,P,N, Sx, bx, Sf, bf, x{k}, xbar, ubar);
    u(k) = un(1);
    x{k+1} = A*x{x} + B*u(k);
    k=k+1;
end

kf = k;


%%

ubar=Nu*rbar % + M*dbar ;
xbar=Nx*rbar

% function un = dlqrcon(A,B,Q,R,P,N,Sx,bx,Su,bu,Sf,bf,x0,xbar,ubar)
un = dlqrcon(A, B, Q, R,P,N, Sx, bx, Sf, bf, x{1}, xbar, ubar);

for k = 1:N
    u(k) = un(k);
    x{k+1} = A*x{k} + B*u(k);
end

for k = N+1:N+30
    u(k) = ubar - K*(x{k} - xbar) 
    x{k+1} = A*x{k} + B*u(k);
end


