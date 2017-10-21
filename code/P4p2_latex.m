%% MATRICES

%System with estimator part IV:
A = [0 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 0;... 
    0 0 0 0 0 1; K_3 0 0 0 0 0];
B = [0 0; 0 K_1; 0 0; K_2 0; 0 0; 0 0];
C = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
D = zeros(3, 2);
sys = ss(A, B, C, D, ...
    'StateName', {'p'; 'p_dot'; 'e'; 'e_dot'; 'lambda'; 'lambda_dot'}, ...
    'InputName', {'V_s'; 'V_d'}, 'OutputName', {'p'; 'e'; 'lambda'});

%Checking observability
O = obsv(sys);
rank(O)


%% LQR from part III

K = lqr(A_L, B_L, Q_L, R_L);
P = inv(C_L*inv(B_L*K-A_L)*B_L);

%% POLE PLACEMENT

system_poles = eig(A_L-B_L*K);

r0 = max(abs(system_poles));

freq = 13;
angle = pi/8;
radius = r0*freq;

spread = -angle:(angle/(2.5)):angle;
observer_poles = -radius * exp(1i*spread);

L = transpose(place(transpose(A),transpose(C),observer_poles));