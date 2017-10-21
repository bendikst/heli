%% PART I
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