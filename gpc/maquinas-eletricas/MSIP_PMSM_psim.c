m = 10^(-3);
k = 10^(3);
fs= 20000;

// modulação
fsw = 20000;
Vdc = 250;
Vin= Vdc;
Vtri = 5;
Vtrip=5;

// Parametros  do motor
P = 6	;	
Polos=P;						
Rs = 0.5;
Ld = 20.1 * m;
Lq = 40.9 * m; 
J = 38.77 * m	;				
Rs= 0.5
Ld= 0.0201
Lq= 0.0409
Vpk_krpm= 278.926066231938

B= 0.0194

Y_Vpk_krpm = 278.926066231938	;
Y = Y_Vpk_krpm * (((2 * pi * k) / 60)^(-1)) ;
lambda_m = ((60)/(sqrt(3)*pi*P*1000))*Y_Vpk_krpm;

//ganhos
kpwm = 1/Vtri;
kr = Vdc;
kiq = 1;
Kni= 1;
Kid= 1;
Kiq= 1;
Kr= 250;

// parametros da malha de controle
Tzid= 1.372043e-4
kcid= 4.369884
Tziq= 1.375231e-4


Tzv= 1.378322e-3
kcv= 1.915463

// Corrente
fc_iq = fsw / 10;
MF_iq_deg = 60;

fc_id = fsw / 10;
MF_id_deg = 60;

// velocidade
fc_mi = fc_iq / (10);
MF_mi_deg = 60;

Tzid= 0.000165152026634352
kcid= 3.62615377522234

Tziq= 0.000165614329878696
kciq= 7.38376564211908

Tzv= 0.0166062920307769
kcv= 0.158983440610366

Vin= 250
Pi= 3.14159265358979


//***********************************************************************************
//***********************************************************************************
// controle de corrente em quadratura (malha interna)
omegaC_iq = 2*pi*fc_iq;
MF_iq_rad = MF_iq_deg * ((pi)/180);

a0_iq = Rs;
a1_iq = Lq;
b0_iq = kpwm*kr*kiq*1;

Modulo_FTLC_iq = (b0_iq)/(sqrt((a0_iq^2)+((omegaC_iq*a1_iq)^2)))
Fase_FTLC_iq = - atan((omegaC_iq*a1_iq)/a0_iq)

//Parametros do controlador PI

omegaZ_iq = omegaC_iq/tan(MF_iq_rad - (pi/2) - Fase_FTLC_iq)
kc_iq = (omegaC_iq/sqrt((omegaC_iq^2)+omegaZ_iq^2))*(1/Modulo_FTLC_iq)

Tz_iq = 1/omegaZ_iq

//***********************************************************************************
//***********************************************************************************
// controle de corrente direta (malha interna)
omegaC_id = 2*pi*fc_id;
MF_id_rad = MF_id_deg * ((pi)/180);

a0_id = Rs;
a1_id = Ld;
b0_id = kpwm*kr*kiq*1;

Modulo_FTLC_id = (b0_id)/(sqrt((a0_id^2)+((omegaC_id*a1_id)^2)))
Fase_FTLC_id = - atan((omegaC_id*a1_id)/a0_id)

//Parametros do controlador PI

omegaZ_id = omegaC_id/tan(MF_id_rad - (pi/2) - Fase_FTLC_id)
kc_id = (omegaC_id/sqrt((omegaC_id^2)+omegaZ_id^2))*(1/Modulo_FTLC_id)

Tz_id = 1/omegaZ_id


//***********************************************************************************
//***********************************************************************************
// controle de velocidade (malha externa)
omegaC_mi = 2*pi*fc_mi;
MF_mi_rad = MF_mi_deg * ((pi)/180);

a1_mi = J;
a0_mi = 0;
b0_mi = (30*3*P*lambda_m)/(pi*2*2);

Modulo_FTLC_mi = (b0_mi)/(sqrt((a0_mi^2)+((omegaC_mi*a1_mi)^2)))
Fase_FTLC_mi = -pi/2

//Parametros do controlador PI

omegaZ_mi = omegaC_mi/tan(MF_mi_rad - (pi/2) - Fase_FTLC_mi)
kc_mi = (omegaC_mi/sqrt((omegaC_mi^2)+omegaZ_mi^2))*(1/Modulo_FTLC_mi)

Tz_mi = 1/omegaZ_mi
