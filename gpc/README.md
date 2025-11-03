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


## GPC

#### 13 May 2025. Generalized Predictive Control for a Single-Phase, Three-Level Voltage Source Inverter Single-phase T-type NPC VSI
Diego Naunay1,Paul Ayala2, Josue Andino3,Wilmar Martinez and Diego Arcos-Aviles2
https://www.mdpi.com/1996-1073/18/10/2541

The proposed control system corresponds to CCS–MPC with a long prediction horizon. Most of the controllers in the literature correspond to FCS–MPC and are based on three-phase grid-connected systems. In contrast to other controllers, the proposed controller is specifically designed for an isolated single-phase VSI and is distinguished by its lower computational cost.
3. Simulation Results
Table 2 presents the simulation parameters used in this study. These parameters are based on a previous study detailed in [42]. Therefore, considering (2), the system’s discretized transfer function, using the ZOH method, at 50 µs is defined in Equation (6), as follows
The controller cost function parameters are tuned offline for simulation purposes, whereas the GPC algorithm is implemented in Matlab® R2023a. This study assumes that the values Np and Nc are equal (Np = Nc = N), the system has no delays (d = 0), and the value δ is set to 1. Since the response with a higher speed and a control action with fewer harmonics is achieved with values of N ≥ 9 [35], a value of N = 9 is considered, as no better performance could be obtained with more extended control and prediction horizons. Additionally, the computational cost of implementing the control algorithm increases as the value of N increases.


#### Basic Problems in Stability and Design of Switched Systems
Liberzon D. and Morse A.S. Basic problems in stability and design of switched systems IEEE Control Syst. Mag. 19 5 59-70 1999
https://experts.illinois.edu/en/publications/basic-problems-in-stability-and-design-of-switched-systems/

#### https://bv.fapesp.br/en/publicacao/75443/control-synthesis-for-dynamic-switched-systems/

#### https://sites.fem.unicamp.br/~grace/ES728_capitulo3.pdf


#### Switched affine systems control design with application to DC–DC converters
Authors: G.S. Deaecto, J.C. Geromel, F.S. Garcia, and J.A. PomilioAuthors Info & Affiliations
Publication: IET Control Theory & Applications
Volume 4, Issue 7
https://doi.org/10.1049/iet-cta.2009.0246


#### Análise de estabilidade e desempenho H2 de sistemas do tipo Lur'e com comutação
DOI: https://doi.org/10.47749/T/UNICAMP.2015.952026


https://sites.fem.unicamp.br/~grace/ES728_capitulo3.pdf
