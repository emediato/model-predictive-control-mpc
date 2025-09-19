function [H11,H22,H12,h,G1,G2] = model_fun(theta1,theta2)

g = 9.8;
L1 = 1;
L2 = 1;
Lc1 = 0.5;
Lc2 = 0.5;
m1 = 20;
m2 = 10;
J1 = 0.8;
J2 = 0.2;

H11 = m1*Lc1^2 + J1 + m2*( L1^2 + Lc2^2 + 2*L1*Lc2*cos(theta2) ) + J2;

H22 = m2*Lc2^2 + J2;

H12 = m2*L1*Lc2*cos(theta2) + m2*Lc2^2 + J2;

h = m2*L1*Lc2*sin(theta2);

G1 = m1*Lc1*g*cos(theta1) + m2*g*( Lc2*cos(theta1 + theta2) + L1*cos(theta1) );

G2 = m2*Lc2*g*cos(theta1 + theta2);
