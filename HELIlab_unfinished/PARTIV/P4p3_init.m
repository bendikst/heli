% FOR HELICOPTER NR 3-10
% This file contains the initialization for the helicopter assignment in
% the course TTK4115. Run this file before you execute QuaRC_ -> Build 
% to build the file heli_q8.mdl.

% Oppdatert høsten 2006 av Jostein Bakkeheim
% Oppdatert høsten 2008 av Arnfinn Aas Eielsen
% Oppdatert høsten 2009 av Jonathan Ronen
% Updated fall 2010, Dominik Breu
% Updated fall 2013, Mark Haring
% Updated spring 2015, Mark Haring


%%%%%%%%%%% Calibration of the encoder and the hardware for the specific
%%%%%%%%%%% helicopter
Joystick_gain_x = 4;
Joystick_gain_y = -7;


%%%%%%%%%%% Physical constants
g = 9.81; % gravitational constant [m/s^2]
l_c = 0.46; % distance elevation axis to counterweight [m]
l_h = 0.66; % distance elevation axis to helicopter head [m]
l_p = 0.175; % distance pitch axis to motor [m]
m_c = 1.92; % Counterweight mass [kg]
m_p = 0.72; % Motor mass [kg]
J_p = 2*m_p * l_p^2;
J_e = m_c * l_c^2 + 2 * m_p * l_h^2;
J_lambda = m_c * l_c^2 + 2 * m_p * (l_h^2 + l_p^2);

K_f = 0.148;
L_1 = l_p * K_f;
L_2 = g*(l_c * m_c - 2 * l_h * m_p);
L_3 = l_h * K_f;
L_4 = -L_3;

K_1 = L_1/J_p;
K_2 = L_3/J_e;
K_3 = -L_2*L_4/(L_3 * J_lambda);


%%%%%%%%%%% Calculated values

K_pp = 14;
K_pd = 2*sqrt(K_pp/K_1);
K_rp = -1.2;

rho = K_rp * K_3;  

%% Matrices
%

A = [0 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 0; 0 0 0 0 0 1; K_3 0 0 0 0 0];
B = [0 0; 0 K_1; 0 0; K_2 0; 0 0; 0 0];
C = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
D = zeros(3, 2);
sys = ss(A, B, C, D, ...
    'StateName', {'p'; 'p_dot'; 'e'; 'e_dot'; 'lambda'; 'lambda_dot'}, ...
    'InputName', {'V_s'; 'V_d'}, 'OutputName', {'p'; 'e'; 'lambda'});

O = obsv(sys);
rank(O)

%% Linear Quadratic Regulator, (From part 3)
%

A_I = [0 1 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 1 0 0 0 0; 0 0 1 0 0];
B_I = [0 0; 0 K_1; K_2 0; 0 0; 0 0];
Q_new = diag([6 3 50 0.01 0.01]); %Changed to reduce integral effect and make power less available
R_I = diag([100, 100]);

K = lqr(A_I, B_I, Q_new, R_I);
%K_P = K(1:2, 1:3);

C_I = [1 0 0 0 0; 0 0 1 0 0]; %Output matrix
P_I = [0 9.999999999999998;9.549869109050652 0];


%% Pole placement
% DESCRIPTIVE TEXT

system_poles = eig(A_I-B_I*K);

r0 = max(abs(system_poles));

freq = 15;
angle = pi/8;
radius = r0*freq;

spread = -angle:(angle/(2.5)):angle;
observer_poles = -radius * exp(1i*spread);

L = transpose(place(transpose(A),transpose(C),observer_poles));

figure(2)
plot(real(system_poles), imag(system_poles), 'sb', real(observer_poles), ... 
    imag(observer_poles), 'rx'); 
grid on; 
axis equal;

%% Additions for problem 4.3
%

%Only measuring e and lambda:
%y = [e; lambda];
C_new = [0 0 1 0 0 0; 0 0 0 0 1 0];
D_new = 0;
sys_new = ss(A, B, C_new, D_new, ...
    'StateName', {'p'; 'p_dot'; 'e'; 'e_dot'; 'lambda'; 'lambda_dot'}, ...
    'InputName', {'V_s'; 'V_d'}, 'OutputName', {'e'; 'lambda'});

new_poles = -[1.5 2 2.5 3 3.5 4];
L_new = transpose(place(transpose(A),transpose(C_new),new_poles));

O_new = obsv(sys_new);
rank(O_new)

%The rank of O_new is also 6.

%Only measuring p and e:
%y = [p, e];

C_bad = [1 0 0 0 0 0; 0 0 1 0 0 0];
sys_bad = ss(A, B, C_bad, D_new, ...
    'StateName', {'p'; 'p_dot'; 'e'; 'e_dot'; 'lambda'; 'lambda_dot'}, ...
    'InputName', {'V_s'; 'V_d'}, 'OutputName', {'e'; 'lambda'});

O_bad = obsv(sys_bad);
rank(O_bad)
%The rank of O_bad is 4, which is not full rank. Hence, the system is not observable!



%% PLOTTING
% Doing some effort to plot this shi*

figure(3)

plot(simout.time, simout.signals.values(:, 13), 'r', simout.time, simout.signals.values(:, 3), 'b--');
legend('Real pitch', 'Estimated pitch')
grid on;


figure(4)

plot(simout.time, simout.signals.values(:, 14), 'r', simout.time, simout.signals.values(:, 4), 'b--');
legend('Real pitch rate', 'Estimated pitch rate')
grid on;


