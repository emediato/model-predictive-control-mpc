sys = ss(-a,b,1,0);
dt = 0.001; nt = (T/dt) + 1;
trec = 0; yrec = 0;
for k = 1:20
    y(k) = yrec(end);
    ...
    t = linspace((k-1)*T,k*T,nt)';
if k < 3
    U = zeros(nt,1);
else
    U = repmat(u(k-2),nt,1); % adding delays
end
yout = lsim(sys,U,t,y(k));
trec = [trec;t(2:end)];
yrec = [yrec;yout(2:end,:)];
end