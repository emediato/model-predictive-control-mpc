H = eye(2);
f = [-5;-5];
b = [2*ones(2,1);0*ones(2,1)];
A = [eye(2);-eye(2)];
% restrição: 0 <= x1 <= 2
bf = quadprog(H,f,A,b)

% B=b;


%% 
a{1} = [0; 0];
Sx = [1 -1 0 0;-1 1 0 0];
Sf = f
N = 15; % Horizonte de controle a ser adotado
%xbar = zeros(4,1); ubar = 0;
xbar = zeros(4,1); ubar = 0;

rbar=1;
un = dlqrcon(A, B, Q, R,P,N, Sx, bx, Sf, bf, a{1}, xbar, ubar);

%%
Gamma = eye(4) ;
max_iter = 10 ;

Sx = [1 -1 0 0;-1 1 0 0];
bx = [0.05;0.05];
bpsi = bx;
Spsi = Sx;


[Sf_1, bf_1, Si_1, bi_1] = determina_oinf(A, Gamma, Spsi, bpsi, max_iter)
