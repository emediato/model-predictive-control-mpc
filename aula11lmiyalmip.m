% addpath(genpath('C:\Users\monic\Desktop\MPClder\matlab\YALMIP-master\YALMIP-master'))
addpath(genpath('C:/Users/monic/Desktop/MPClder/matlab/YALMIP-master/YALMIP-master'))
% LINEAR MATRIX INEQUALITY
% LMI
% P > 0
% (A - BK)_TRANSPOSE . P(A - BK) - P < 0 
sdpvar x % Vari´avel escalar
restricoes = [ x <= 5 ];
custo = (x - 10)^2;
optimize(restricoes,custo);
solucao = value(x)

A = [1 1;0 2]; B = [0;1]; K = [0 0];
eig(A)
% 2 >1 > INSTABLE
%% 
n = size(A,1);
P = sdpvar(n); % Vari´avel matricial (n x n) sim´etrica
restricoes = [ P >= 0 , (A-B*K)'*P*(A-B*K)-P <= 0 ];
ops = sdpsettings('solver','sedumi')
optimize(restricoes,[],ops)
%%
Psol = value(P)
eig(Psol)
eig((A-B*K)'*Psol*(A-B*K) - Psol)