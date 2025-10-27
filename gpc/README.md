## Models from 

rojas2021.pdf
chirantan2020-2.pdf

## Support

https://github.com/treezao/Livro_mpc_codes/tree/main/Volume1/Capitulo4-GPC

## space vector modulation

Vdc = ...;        % DC link voltage
d = ...;          % Modulation index
wo = ...;         % Output angular frequency
t = ...;          % Time vector

va = d * Vdc * cos(wo * t);
vb = d * Vdc * cos(wo * t - 2*pi/3);
v_cz_star = d * Vdc * cos(wo * t + 2*pi/3);

v_alpha = (2/3) * va - (1/3) * vb - (1/3) * vc;
v_beta  = (sqrt(3)/3) * (vb - vc)
V_ref = sqrt(v_alpha^2 + v_beta^2);
theta = atan2(v_beta, v_alpha);



