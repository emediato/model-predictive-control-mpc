%aula9

n = 4; q = 1; p = 1;
Achif = [Af B*(K*Nx + Nu);zeros(q,n) eye(q)];
Gamma = [eye(n) zeros(n,q);
-K (K*Nx + Nu);
zeros(n,n) Nx;
zeros(p,n) Nu];
Spsi = blkdiag(Sx,Su,Sx,Su);


% blkdiag for diag, horzcat, vertcat
% spdiags, triu, tril, blkdiag.

epsilonx = 1e-3*ones(size(Sx,1),1);
epsilonu = 1e-3*ones(size(Su,1),1);
bpsi = [bx;bu;bx-epsilonx;bu-epsilonu];

 max_iter = 100; tol = 1e-6;
[So,bo] = determina_oinf(Achif,Gamma,Spsi,bpsi,max_iter,tol);


rbar = pi;
Ox = So(:,1:n);
Or = So(:,n+1:end);
Sf = Ox;
bf = bo - Or*rbar;
Oinf = Polyhedron('H',[Sf bf]);
plot(Oinf)