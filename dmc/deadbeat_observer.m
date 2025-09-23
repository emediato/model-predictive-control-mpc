E = B;
Achi = [A E;zeros(2,4) eye(2)];
Bchi = [B;zeros(2,2)];
H = eye(4);
Hchi = [H zeros(4,2)];

eig_des = [0:5]*1e-2;

L = (place(Achi’,Hchi’,eig_des))’;
M = inv(Achi)*L;
eig(Achi - M*Hchi*Achi)
aux = inv([A - eye(4) B;C zeros(2,2)])*[-E;zeros(2,2)];
Mx = aux(1:4,:)
Mu = aux(5:6,:);
xnl{1} = xnlbar;
trec = 0; xnlrec = xnl{1}’;
chi_pri{1} = zeros(6,1);


for k = 1:200
  x{k} = xnl{k} - xnlbar;
  zpri{k} = Hchi*chi_pri{k};
  z{k} = x{k};
  chi_pos{k} = chi_pri{k} + M*(z{k} - zpri{k});
  dbar_pos{k} = chi_pos{k}(5:6);
  xbar{k} = Nx*rbar + Mx*dbar_pos{k};
  ubar{k} = Nu*rbar + Mu*dbar_pos{k};
  u{k} = ubar{k} - K*(x{k} - xbar{k});
  chi_pri{k+1} = Achi*chi_pos{k} + Bchi*u{k};
  unl{k} = u{k} + unlbar;
  [tout,xnlout] = ode45(@(t,xnl) ...
  equacao_estado(t,xnl,unl{k}),[(k-1)*T k*T],xnl{k});
  trec = [trec;tout(2:end)];
  xnlrec = [xnlrec;xnlout(2:end,:)];
  xnl{k+1} = xnlrec(end,:)’;
end
