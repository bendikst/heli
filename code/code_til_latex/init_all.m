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


%% PHYSICAL CONSTANTS
%Defining given constants and the constants we have derived.

g = 9.81; % gravitational constant [m/s^2]
l_c = 0.46; % distance elevation axis to counterweight [m]
l_h = 0.66; % distance elevation axis to helicopter head [m]
l_p = 0.175; % distance pitch axis to motor [m]
m_c = 1.92; % Counterweight mass [kg]
m_p = 0.72; % Motor mass [kg]
J_p = 2*m_p * l_p^2; %Moment of inertia
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


%% CALCULATED VALUES
%Part II


K_f = 0.148;
L_1 = l_p * K_f;
L_2 = g*(l_c * m_c - 2 * l_h * m_p);
L_3 = l_h * K_f;
L_4 = -L_3;

K_1 = L_1/J_p;
K_2 = L_3/J_e;
K_3 = -L_2*L_4/(L_3 * J_lambda);

K_pp = 14;
K_pd = 2*sqrt(K_pp/K_1);
K_rp = -1.2;

rho = K_rp * K_3;  

%% Matrices
%Part III, problem 2
A_L = [0 1 0; 0 0 0; 0 0 0];
B_L = [0 0; 0 K_1; K_2 0];
C_L = [1 0 0; 0 0 1];
Q_L = diag([91.2 50 100]);
R_L = diag([1 1]);

%Part III, problem 3
A_I = [0 1 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 1 0 0 0 0; 0 0 1 0 0];
B_I = [0 0; 0 K_1; K_2 0; 0 0; 0 0];
C_I = [1 0 0 0 0; 0 0 1 0 0]; %Output matrix
Q_I = diag([91.2 50 100 50 50]);
R_I = diag([1, 1]);

%Part IV
A_E = [0 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 0;...
    0 0 0 0 0 1; K_3 0 0 0 0 0];
B_E = [0 0; 0 K_1; 0 0; K_2 0; 0 0; 0 0];
C_E = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
D_E = zeros(3, 2);
sys = ss(A_E, B_E, C_E, D_E, ...
    'StateName', {'p'; 'p_dot'; 'e'; 'e_dot'; 'lambda'; 'lambda_dot'}, ...
    'InputName', {'V_s'; 'V_d'}, 'OutputName', {'p'; 'e'; 'lambda'});

O = obsv(sys);
rank(O)

%% Linear Quadratic Regulator
%Part III

K_L = lqr(A_L, B_L, Q_L, R_L);
P_L = inv(C_L*inv(B_L*K_L-A_L)*B_L);

K_I = lqr(A_I, B_I, Q_I, R_I);
P_I = [0 9.999999999999998;9.549869109050652 0];


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
