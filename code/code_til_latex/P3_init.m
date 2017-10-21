%% Matrices
%Part III - problem 2
A_L = [0 1 0; 0 0 0; 0 0 0];
B_L = [0 0; 0 K_1; K_2 0];
C_L = [1 0 0; 0 0 1];
Q_L = diag([91.2 50 100]);
R_L = diag([1 1]);

%Part III - problem 3 - integrator
A_I = [0 1 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 1 0 0 0 0; 0 0 1 0 0];
B_I = [0 0; 0 K_1; K_2 0; 0 0; 0 0];
C_I = [1 0 0 0 0; 0 0 1 0 0];
Q_I = diag([91.2 50 100 50 50]);
R_I = diag([1, 1]);

%% Linear Quadratic Regulator
%Part III - problem 2

K_L = lqr(A_L, B_L, Q_L, R_L);
P_L = inv(C_L*inv(B_L*K_L-A_L)*B_L);

%Part III - problem 3
K_I = lqr(A_I, B_I, Q_I, R_I);
P_I = [0 9.999999999999998;9.549869109050652 0];