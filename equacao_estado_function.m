function xnldot = equacao_estado_function(t,xnl,unl);

theta1 = xnl(1); theta2 = xnl(2);
theta1dot = xnl(3); theta2dot = xnl(4);

tau1 = unl(1); tau2 = unl(2);
[H11,H22,H12,h,G1,G2] = model_fun(theta1,theta2);
H = [H11 H12;H12 H22];


aux = inv(H)*[tau1 + h*(theta2dot)^2 + 2*h*theta1dot*theta2dot - G1;
              tau2 - h*(theta1dot)^2 - G2];

xnldot = zeros(4,1);          
xnldot(1) = xnl(3);
xnldot(2) = xnl(4);
xnldot(3) = aux(1);
xnldot(4) = aux(2);