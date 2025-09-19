
function un = dlqrcon(A,B,Q,R,P,N,Sx,bx,Su,bu,Sf,bf,x0,xbar,ubar)

Rn = R; Sun = Su;
An = A;
for i = 2:N
    Rn = blkdiag(Rn,R);
    Sun = blkdiag(Sun,Su);
    An = [An;A^i];
end
size(An)
bun = repmat(bu,N,1);
n = size(B,1);
p = size(B,2);

%  mat2cell Break matrix up into a cell array of matrices.
%     C = mat2cell(X,M,N) breaks up the 2-D array X into a cell array of  
%     adjacent submatrices of X. X is an array of size [ROW COL], M is the 
%     vector of row sizes (must sum to ROW) and N is the vector of column 
%     sizes (must sum to COL). The elements of M and N determine the size of
%     each cell in C by satisfying the following formula for I = 1:LENGTH(M)
%     and J = 1:LENGTH(N),

Baux = mat2cell(zeros(n*N,p*N),n*ones(N,1),p*ones(1,N));
for i = 1:N
    for j = 1:i
        Baux{i,j} = (A^(i-j))*B;
    end
end
Bn = cell2mat(Baux);

Qn = Q; Sxn = Sx;
for i = 2:N-1
    Qn = blkdiag(Qn,Q);
    Sxn = blkdiag(Sxn,Sx);
end
Qn = blkdiag(Qn,P);
Sxn = blkdiag(Sxn,Sf);
bxn = [repmat(bx,N-1,1);bf];
size(bxn)
size(Sxn)
size(bun)
size(Sxn*An*x0)
size(bxn - Sxn*An*x0)

% Uso do quadprog
Aqp = [Sxn*Bn;Sun];
bqp = [bxn - Sxn*An*x0;bun];

Hqp = Bn'*Qn*Bn + Rn;

% Simetrização
Hqp = (Hqp + Hqp')/2;

% max(max(abs(Hqp - Hqp')))

fqp = (Bn'*Qn*(An*x0 - repmat(xbar,N,1)) - Rn*repmat(ubar,N,1));

% options = optimset('Algorithm','active-set');
% un = quadprog(Hqp,fqp,Aqp,bqp,[],[],[],[],[],options);

un = quadprog(Hqp,fqp,Aqp,bqp);
