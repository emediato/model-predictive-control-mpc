ppppppppppppppA = [1 0.1;0 1]; B = [5e-3;0.1];
Q = diag([50,1]); R = 1;
K = dlqr(A,B,Q,R);
C = [1 0];

aux = inv([A - eye(2) B;C 0])*[0;0;1]; % nao considera perturbação
Nx = aux(1:2)
Nu = aux(3)

rbar = pi;
xbar = Nx*rbar; ubar = Nu*rbar;

Gamma = eye(4) ;
max_iter = 100 ;
tol = 1e-6;

Sx = [1 -1 0 0;-1 1 0 0];
bx = [0.05;0.05];
bpsi = bx;
Spsi = Sx;

Af = A - B*K

n = 2; q = 1; p = 1;
Achif = [Af B*(K*Nx + Nu);
  
[Sf, bf, Si, bi] = determina_oinf(Achif,Gamma,Spsi,bpsi,max_iter);

%%
Ox = So(:,1:n); Or = So(:,n+1:end);
Sf = Ox; bf = bo - Or*rbar;

n = size(B,1); p = size(B,2); q = size(So,2) - n;

x{1} = [0;0];
for k = 1:40
    u(k) = ubar - K*(x{k} - xbar);
    x{k+1} = A*x{k} + B*u(k);
end

[Sf, bf, Si, bi] = determina_oinf(A, Gamma, Spsi, bpsi, max_iter)


unk = dlqrcon(A,B,Q,R,P,N,Sx,bx,Su,bu,Sf,bf,xk,xbar,ubar);

p = size(B,2);
sys = unk(1:p); % Sa´ıda da S-function: u(kT)


