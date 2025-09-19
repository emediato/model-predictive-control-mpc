Ac=[0 1; 852 0]
Bc=[0; -20.5]
T=0.005;

m1 = 5; m2 = 3; ks = 10; g = 9.8;
Ac = [0 0 1 0;
      0 0 0 1;
      -ks/m1 ks/m1 0 0;
      ks/m2 -ks/m2 0 0];
Bc = [0;0;1/m1;0];
Ec = [0;0;-1;-1];
C = [0 1 0 0];


[A,B]= c2dm(Ac, Bc, [], [], T, 'zoh')

C=[1 0]
T = 0.2; % sample
%angulo inclinação
x2=0.05;

dbar = g * sin(x2) % perturbação cte


% E=Ec * dbar %%%%%%%%%%%%%%%%%%%%%%%%% verify
Sysc = ss(Ac, Bc, [], []);
Sys=c2d(Sysc, T, 'zoh');

A=Sys.A
B=Sys.B

[~, E] =c2dm(Ac, Ec, [], [], T, 'zoh')


% autovalores sao as raizes de um polinomio caracteristico
xi=0.69; % considerando dinamica de freq natural
wn=1.8/6 ; % 6 segundos 
% controle nao deve ultrapassar 2N

eigc=roots([1 2*xi*wn wn^2])

eigd=exp(eigc*T)
J=eig(A-B*K(1))

K=acker(A,B,J) 
k1=place(A,B,J) 


aux = inv([A - eye(4) B; C 0])*[0; 0; 0; 0; 1]

Nx = aux(1:4)
Nu = aux(5)

rbar=1e-1;
ubar=Nu*rbar % + M*dbar ;
xbar=Nx*rbar % + Mx*dbar ;

u=[0;0;0;0];
trec = 0;

nt=30
x= repmat([0;0;0;0]',nt,1)
u = [1;1;1;1];
U = repmat(u',nt,1)


for k = 1:30
    x{k} = x{k+1} - xnlbar;
    
    u{k} = ubar -K*(x{k} - xbar);
    
    unl{k} = u{k} + unlbar;
    
    [tout,  xnlout] 
    trec = [trec;  tout(2:end)];
    xnlrec = [xnlrec;xnlout(2:end,:)];
    
    x{k+1} = A*x{k} + B*u{k};

    %xnl{k+1} = xnlrec(end,:)';
end



