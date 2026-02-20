beta_tt = 2;
m_dot = 3;
Pt1 = 1e5;
Tt1 = 300;
cp_exp = 0.037; cv_exp = 0.028; Mm = 17.03;
cp = cp_exp/Mm * 1e6; cv = cv_exp/Mm * 1e6;
gamma = cp_exp/cv_exp; R = cp-cv;
omega = 15000 * 2*pi/60;

dht_is = gamma * R/(gamma -1) * Tt1 * (beta_tt ^((gamma - 1)/gamma) -1);
rho_t1 = Pt1/R/Tt1;
Q_in = m_dot/rho_t1;
Ds = 5.0;
D2 = Ds * sqrt(Q_in)/(dht_is^(1/4));
L = gamma * R/(gamma -1) * Tt1/0.82 * (beta_tt ^((gamma - 1)/gamma) -1);
fprintf('L = %.2f, D2 = %.4f
', L, D2);
