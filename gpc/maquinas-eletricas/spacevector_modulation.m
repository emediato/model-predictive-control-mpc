Vdc = 400;
switch_state = [
    0 0 0;
    0 0 1;
    0 1 0;
    0 1 1;
    1 0 0;
    1 1 0;
    1 0 1;
    1 1 1;
];

v_alpha = zeros(size(switch_state,1),1);
v_beta = zeros(size(switch_state,1),1);

for k = 1:size(switch_state,1)
    Sa = switch_state(k,1);
    Sb = switch_state(k,2);
    Sc = switch_state(k,3);
    
    va = (2*Sa - 1)*Vdc/2;
    vb = (2*Sb - 1)*Vdc/2;
    vc = (2*Sc - 1)*Vdc/2;
    
    v_alpha(k) = (2/3)*(va - 0.5*vb - 0.5*vc);
    v_beta(k)  = (2/3)*((sqrt(3)/2)*(vb - vc));
end

% Compute modulation index for each vector
Vref = sqrt(v_alpha.^2 + v_beta.^2);
Vmax = Vdc/sqrt(3);
m = Vref / Vmax;

disp('Modulation index for each switching state:');
disp(m);

% Optional: plot space vectors
figure;
plot(v_alpha, v_beta, 'o');
xlabel('v_\alpha'); ylabel('v_\beta');
title('Space Vectors from Switch States');
grid on;
