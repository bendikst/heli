%% Matrices
%Part IV

A_E = [0 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 0;...
    0 0 0 0 0 1; K_3 0 0 0 0 0];
B_E = [0 0; 0 K_1; 0 0; K_2 0; 0 0; 0 0];
C_E = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
D_E = zeros(3, 2);
sys = ss(A_E, B_E, C_E, D_E, ...
    'StateName', {'p'; 'p_dot'; 'e'; 'e_dot'; 'lambda'; 'lambda_dot'}, ...
    'InputName', {'V_s'; 'V_d'}, 'OutputName', {'p'; 'e'; 'lambda'});

%Checking observability
O = obsv(sys);
rank(O)

%% Pole placement
% Part IV - problem 2

system_poles_L = eig(A_L-B_L*K_L);
system_poles_I = eig(A_I-B_I*K_I);

r0_L = max(abs(system_poles_L));
r0_I = max(abs(system_poles_I));

freq = 13;
angle = pi/8;
spread = -angle:(angle/(2.5)):angle;

radius_L = r0_L*freq;
radius_I = r0_I*freq;

observer_poles_L = -radius_L * exp(1i*spread);
observer_poles_I = -radius_I * exp(1i*spread);

L_L = transpose(place(transpose(A_E),transpose(C_E),observer_poles_L));
L_I = transpose(place(transpose(A_E),transpose(C_E),observer_poles_I));

%% Additions for problem 4.3
%
%Only measuring e and lambda:

C_new = [0 0 1 0 0 0; 0 0 0 0 1 0];
D_new = 0;
sys_new = ss(A_E, B_E, C_new, D_new, ...
    'StateName', {'p'; 'p_dot'; 'e'; 'e_dot'; 'lambda'; 'lambda_dot'}, ...
    'InputName', {'V_s'; 'V_d'}, 'OutputName', {'e'; 'lambda'});

new_poles = -[1.5 2 2.5 3 3.5 4];
L_new = transpose(place(transpose(A_E),transpose(C_new),new_poles));

O_new = obsv(sys_new);
rank(O_new)

%The rank of O_new is also 6.