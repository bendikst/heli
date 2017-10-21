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