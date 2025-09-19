ops = sdpsettings('solver','sedumi');
optimize(lmis,g,ops)
Ysol = value(Y);
Xsol = value(X);
K = -Ysol*inv(Xsol);