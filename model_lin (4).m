function [Ac,Bc] = model_lin(theta1,theta2)

g = 9.8;
L1 = 1;
Lc1 = 0.5;
Lc2 = 0.5;
m1 = 20;
m2 = 10;
J1 = 0.8;
J2 = 0.2;

a1 = m1*Lc1^2 + J1 + m2*(L1^2 + Lc2^2 + 2*L1*Lc2*cos(theta2)) + J2;
a2 = m2*L1*Lc2*cos(theta2) + m2*Lc2^2 + J2;
a3 = (m1*Lc1 + m2*L1)*g*sin(theta1) + m2*Lc2*g*sin(theta1 + theta2);
a4 = m2*Lc2*g*sin(theta1 + theta2);
a5 = m2*L1*Lc2*cos(theta2) + m2*Lc2^2 + J2;
a6 =  m2*Lc2^2 + J2;
a7 = m2*Lc2*g*sin(theta1 + theta2);
a8 = a7;

aux1 = inv([a1 a2;a5 a6]);
aux2 = aux1*[a3 a4;a7 a8];

Ac = [0 0 1 0;
      0 0 0 1;
      aux2 zeros(2,2)];

Bc = [zeros(2,2);aux1];