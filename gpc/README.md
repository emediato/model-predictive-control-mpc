## Models from 

rojas2021.pdf
chirantan2020-2.pdf

## Support

https://github.com/treezao/Livro_mpc_codes/tree/main/Volume1/Capitulo4-GPC

## PWM VSI modulation
https://www.youtube.com/watch?v=rIer4hAUBfk

## space vector modulation
https://www.youtube.com/watch?v=kY6v8PIFePo

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

https://www.researchgate.net/publication/346027016_DISSERTACAO_-_Analise_e_Concepcao_de_um_Sistema_Hibrido_de_Armazenamento_de_Energia_para_Aplicacao_em_Locomotivas_Diesel-Eletricas


## GPC

#### 13 May 2025. Generalized Predictive Control for a Single-Phase, Three-Level Voltage Source Inverter Single-phase T-type NPC VSI
Diego Naunay1,Paul Ayala2, Josue Andino3,Wilmar Martinez and Diego Arcos-Aviles2
https://www.mdpi.com/1996-1073/18/10/2541

#### Generalised predictive controller (GPC) design on single-phase full-bridge inverter with a novel identification method
Mehran Jelodari Memeghani, Seyyed Morteza Ghamari, Taha Yousefi Jouybari, Hasan Mollaee, Patrick Wheeler
First published: 23 May 2022 https://doi.org/10.1049/cth2.12295
https://ietresearch.onlinelibrary.wiley.com/doi/full/10.1049/cth2.12295
The proposed control system corresponds to CCS–MPC with a long prediction horizon. Most of the controllers in the literature correspond to FCS–MPC and are based on three-phase grid-connected systems. In contrast to other controllers, the proposed controller is specifically designed for an isolated single-phase VSI and is distinguished by its lower computational cost.
Simulation Results
Table 2 presents the simulation parameters used in this study. These parameters are based on a previous study detailed in [42]. Therefore, considering (2), the system’s discretized transfer function, using the ZOH method, at 50 µs is defined in Equation (6), as follows
The controller cost function parameters are tuned offline for simulation purposes, whereas the GPC algorithm is implemented in Matlab® R2023a. This study assumes that the values Np and Nc are equal (Np = Nc = N), the system has no delays (d = 0), and the value δ is set to 1. Since the response with a higher speed and a control action with fewer harmonics is achieved with values of N ≥ 9 [35], a value of N = 9 is considered, as no better performance could be obtained with more extended control and prediction horizons. Additionally, the computational cost of implementing the control algorithm increases as the value of N increases.


The GMPC method is a well-known MPC method since it contains an excellent performance with a great robustness for both academia and industry. The principle behind GMPC is based on calculating a signal sequence of future control, while the goal is to reach a minimum value of a multistage cost function depicted on a prediction horizon [17]. In addition, the optimisation index is the expectation of a quadratic functions calculating the difference between the predicted reference and the predicted system outputs; also, this reference is obtained on the horizon, where the control effort is measured by a quadratic function. The general block diagram of GMPC is shown in Figure 3.
The design of a GPC on single-phase full-bridge inverter is presented by this work using an LC filter to control the current feed to the load. The efficiency of the controller has been tested by simulation and experimental results. This predictive algorithm is classified as CCS-MPC, which needs an accurate parametric estimation. 
EXPERIMENTAL AND SIMULATION RESULTS
It should be noted that zero-order hold (ZOH) technique is utilised to discrete the controller. Also, the PWM scheme is used to fire the converter switches. To clarify this technique for better understanding, one can assume that the PWM is a block that gets the control signal and compare it with a triangular signal, then a square wave is generated that has the role of actuating the switches. In Figure 5a hardware real-time implementation of converter is shown consisting of the power and control realisation, simultaneously.


## FCS

#### A Very Simple Strategy for High Quality Performance of AC Machines Using Model Predictive Control
March 2018IEEE Transactions on Power Electronics PP(99):1-1
DOI:10.1109/TPEL.2018.2812833



## Control

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
