
clear;
clc;
m = 10^(-3);
k = 10^(3);

// Parametros  do motor
P = 6	;							
Rs = 0.5;
Ld = 20.1 * m;
Lq = 40.9 * m; 
J = 0.03877	;				
Y = 0.5126 ;
Y_Vpk_krpm = 53.679346	;

lamnda_m = ((60)/(sqrt(3)*pi*P*1000))*Y_Vpk_krpm;

// modulação
fsw = 20000;
Vdc = 96;
Vtri = 5;

kpwm = 1/Vtri;
kr = Vdc;
kiq = 1;

// parametros da malha de controle
// Corrente
fc_iq = fsw / 10;
MF_iq_deg = 60;

// velocidade
fc_mi = fc_iq / 50;
MF_mi_deg = 60;


//*************************************************//
// controle de corrente (malha interna)
wc_iq = 2*pi*fc_iq;
MF_iq_rad = MF_iq_deg * ((pi)/180);

a0_iq = Rs;
a1_iq = Lq;
b0_iq = kpwm*kr*kiq*1;

Modulo_FTLC_iq = (b0_iq)/(sqrt((a0_iq^2)+((wc_iq*a1_iq)^2)))
Fase_FTLC_iq = - atan((wc_iq*a1_iq)/a0_iq)

//Parametros do controlador PI

wz_iq = wc_iq/tan(MF_iq_rad - (pi/2) - Fase_FTLC_iq)
kc_iq = (wc_iq/sqrt((wc_iq^2)+wz_iq^2))*(1/Modulo_FTLC_iq)

Tz_iq = 1/wz_iq

//*************************************************//
// controle de velocidade (malha externa)

wc_mi = 2*pi*fc_mi;
MF_mi_rad = MF_mi_deg * ((pi)/180);


a1_mi = J;
a0_mi = 0;
b0_mi = (30*3*P*lamnda_m)/(pi*2*2);

Modulo_FTLC_mi = (b0_mi)/(sqrt((a0_mi^2)+((wc_mi*a1_mi)^2)))
Fase_FTLC_mi = -pi/2

//Parametros do controlador PI

wz_mi = wc_mi/tan(MF_mi_rad - (pi/2) - Fase_FTLC_mi)
kc_mi = (wc_mi/sqrt((wc_mi^2)+wz_mi^2))*(1/Modulo_FTLC_mi)

Tz_mi = 1/wz_mi
