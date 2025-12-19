clear
close all
clc
%% input
beta_tt = 2;
m_dot = 3; %kg/s
Pt1 = 1 * 1e5; % Pa
Tt1 = 300; %K
%%% properties gas %%%
cp_exp = 0.037; 
cv_exp = 0.028; 
Mm = 17.03; 
cp = cp_exp/Mm * 1e6;
cv = cv_exp/Mm * 1e6;
gamma = cp_exp/cv_exp;
R = cp-cv; 
%%% Sutherland's Law for Ammonia)
mu_0_NH3 = 9.82e-6;   
T_0_NH3  = 293.15;    
S_NH3    = 370;       
mu_NH3 = @(T) mu_0_NH3 * (T./T_0_NH3).^(1.5) .* (T_0_NH3 + S_NH3) ./ (T + S_NH3);
%%% assumptions %%%
rpm = 15000; 
omega = rpm * 2*pi/60;
alpha_2 = 70 * pi/180;
%% balje
dht_is = gamma * R/(gamma -1) * Tt1 * (beta_tt ^((gamma - 1)/gamma) -1);
rho_t1 = Pt1/R/Tt1;
Q_in = m_dot/rho_t1;
omega_s = omega * sqrt(Q_in)/dht_is^(3/4);
Ds = 4.5;  
D2 = Ds * sqrt(Q_in)/(dht_is^(1/4));
rho_vect = [rho_t1];
err_new = 1;
iter = 0;
relaxation = 0.8;
while abs(err_new) > 1e-5
    D1_h = 0.2; 
    R1_h = D1_h/2;
    rho = rho_vect(end);
    cm1_max = sqrt(2*cp*Tt1);
    r_min = sqrt(m_dot/rho/cm1_max/pi + R1_h^2);
    R1_t = linspace(r_min*1.000001,R1_h*8, 100000);
    cm1 = @(r) m_dot ./ (rho * pi * (r.^2 - R1_h^2));
    T1 = @(r) Tt1 - cm1(r).^2/(2*cp);
    M1_tip= @(r) sqrt(((omega.*r).^2 + cm1(r).^2) ./ (gamma * R * T1(r)));
    M1_t_eval = arrayfun(M1_tip,R1_t);
    [M1_tip_rel, idx] = min(M1_t_eval);
    R1_t_opt = R1_t(idx);
    T1 = Tt1 - cm1(R1_t_opt).^2/(2*cp);
    P1 = Pt1 * (T1/Tt1)^(gamma/(gamma-1));
    rho_new = P1/R/T1;
    rho_new = rho_vect(end) + relaxation * (rho_new - rho_vect(end));
    err_new = abs(rho_new - rho)/rho;
    rho_vect = [rho_vect; rho_new];
    iter = iter+1;
end
rho1 = rho_vect(end);

V1_ref = m_dot / (rho1 * pi * (R1_t_opt^2 - R1_h^2));
R1_m = (R1_t_opt + R1_h)/2;
D1_m = R1_m * 2;
b1 = R1_t_opt - R1_h;
A1 = pi * D1_m * b1;
D1_t_opt = R1_t_opt * 2;
U1_mean = omega * R1_m;
U1_tip = omega * R1_t_opt;
U1_hub = omega * R1_h;
W1_mean = sqrt(U1_mean^2 + V1_ref^2);
W1_meridional = V1_ref; 
W1_tip_t = -U1_tip; W1_hub_t = -U1_hub;
W1_tip = sqrt(W1_tip_t^2 + W1_meridional^2);
W1_hub = sqrt(W1_hub_t^2 + W1_meridional^2);
beta1_mean = atan(-U1_mean / V1_ref);
beta1_tip = atan(W1_tip_t/W1_meridional);
beta1_hub = atan(W1_hub_t/W1_meridional);
V1 = V1_ref; 
T1 = Tt1 - V1^2/(2*cp);

%% SETUP SENSITIVITY
n_sens = 50;
D2_vec = linspace(D2 * 0.8, D2* 1.2, n_sens);
eta_sens = zeros(1, n_sens);
b2_sens  = zeros(1, n_sens);

for i = 1:n_sens
    D2_curr = D2_vec(i);
    R2_curr = D2_curr / 2;
    U2_curr = omega * R2_curr;
    eta_tt_guess = 0.85;
    rho2_guess = rho1 * 1.5; 
    err_loc = 1;
    iter_loc = 0;
    
    while abs(err_loc) > 1e-5 && iter_loc < 100
        iter_loc = iter_loc + 1;
        L_req = cp * Tt1/eta_tt_guess * (beta_tt^((gamma-1)/gamma) - 1);
        V2_tg_curr = L_req / U2_curr;
        V2_m_curr = V2_tg_curr / tan(alpha_2);
        V2_curr = sqrt(V2_m_curr^2 + V2_tg_curr^2);
        W2_tg_curr = V2_tg_curr - U2_curr;
        W2_curr = sqrt(V2_m_curr^2 + W2_tg_curr^2);
        Tt2_curr = Tt1 + L_req/cp;
        T2_curr = Tt2_curr - V2_curr^2/(2*cp);
        M2_curr = V2_curr / sqrt(gamma * R * T2_curr);
        Pt2_curr = Pt1 * (1 + eta_tt_guess * L_req/(cp*Tt1))^(gamma/(gamma-1));
        P2_curr = Pt2_curr / (1 + (gamma-1)/2 * M2_curr^2)^(gamma/(gamma-1));
        rho2_new = P2_curr / (R * T2_curr);
        rho2_guess = 0.5 * rho2_guess + 0.5 * rho2_new;
        b2_curr = m_dot / (rho2_guess * V2_m_curr * pi * D2_curr);
        beta2_curr = atan(W2_tg_curr / V2_m_curr);
        beta_av_curr = 0.5 * (beta1_mean + beta2_curr);
        N_bl_curr = ceil(2*pi*cos(beta_av_curr)/(0.4*log(D2_curr/(D1_m))));
        if mod(N_bl_curr,2)==1
            N_bl_curr = N_bl_curr + 1;
        end
        t_curr = max(0.003 * D2_curr, 0.0025);
        
        D2 = D2_curr; R2 = R2_curr; U2 = U2_curr;
        b2 = b2_curr; N_bl = N_bl_curr; t = t_curr;
        V2_tg = V2_tg_curr; V2_meridional = V2_m_curr; V2 = V2_curr;
        W2 = W2_curr; beta2 = beta2_curr;
        T2 = T2_curr; L = L_req; rho2 = rho2_guess;
        L_eul = U2 * V2_tg;
        phi = 4 * Q_in/(pi*U2*D2^2);
        beta1_geom_mean = atan(1 - t*N_bl/(2*pi*R1_m)*tan(beta1_mean));
        beta1_geom_tip = atan(1 - t*N_bl/(2*pi*R1_t_opt)*tan(beta1_tip)); 
        mu_slip = 1 - 0.63*pi/N_bl;
        V2_tg_inf = (1- mu_slip)*U2 + V2_tg;
        beta2_geom = atan((V2_tg_inf - U2)/V2_meridional);
        solidity = 1/0.4;
        dtheta = log(D2/D1_m) * tan(beta2);
        c = R2 * dtheta/(sin(beta2));
        s = c/solidity;
        
        % Mapping variabili inlet mancanti per il blocco perdite
        W1_mean_t = -U1_mean; 
        
        dH_inc = 0.4*(W1_mean - V1/cos(beta1_mean))^2;
        %%% IMPELLER INTERNAL
        % skin friction (Jansen, 1967)
        visc_din1 = mu_NH3(T1); %sutherland
        % D_h = 4*(2*pi*R2*b2/N_bl)/(2*b2+2*2*pi*R2/N_bl);
        D_h = D2*cos(beta2) / (( N_bl/pi + D2*cos(beta2)/b2 ) + ( 0.5*(D1_t_opt/D2 + D1_h/D2) * ((cos(beta1_tip) + cos(beta1_hub))/2) )/( N_bl/pi + ((D1_t_opt + D1_h)/(D1_t_opt - D1_h)) * ((cos(beta1_tip) + cos(beta1_hub))/2) ) );
        Re_f = U2 * D_h/visc_din1*rho_t1;
        Cf = 0.0412 * Re_f ^(-0.1925);
        Lz = 0.08 + 3.16*phi; %computed as aungier
        %Lz = 0.1 + 2*phi by Hamid Hazby (sembra meglio sperimentalmente)
        lb = pi/8 *(D2- 1/2 * (D1_t_opt+D1_h) - b2 + 2*Lz)*(2/((cos(beta1_tip)+cos(beta1_hub))/2 + (cos(beta2))));
        W12_av = 1/4 * (W1_hub + W1_tip + 2*W2);
        dH_fr = (2 * Cf * W12_av^2 * lb)/(D_h);
        % clearance (Jansen)
        eps = 5e-4;
        dH_cl = 0.6 * eps / b2 * V2_tg * ( 4*pi / (b2 * N_bl) * (R1_t_opt^2 - R1_h^2) / (R2 - R1_t_opt) *V2_tg * V1 / (1 + rho2 / rho1) )^(1/2);
        % blade loading (Aungier, 1995) assumption of no inlet swirl
        dW = 2*pi*D2*V2_tg/(N_bl*lb);
        dH_bl = 1/48*dW^2;
        % wake mixing (Aungier, 1995)
        W_max = 0.5*(W1_mean + W2 + dW);
        D_eq = W_max/W2;
        if D_eq <= 2
            W_sep = W2;
        else
            W_sep = W2*D_eq*0.5;
        end
        t_te = 2e-3; % imposto da geometria
        eps2 = 1- (N_bl * t_te)/(pi*D2);
        A2 = pi*b2*D2*eps2;
        W_out = sqrt((V2_meridional*A2/(pi*D2*b2))^2+W1_mean_t^2);
        dH_mix = 1/2 * (W_sep - W_out)^2;
        % entrance diffusion (Aungier)
        % pitch = pi*D1_m/N_bl - t;
        A_th_geom = s * cos(beta1_geom_mean)*b1;
        A_th = 0.97 * A_th_geom;
        W_th = m_dot/N_bl/rho1/A_th;
        dH_diff = 0.4*(W1_mean - W_th)^2;
        if (W1_tip/W_th > 1.75) && (dH_diff < 0.5*(W1_tip - 1.75*W_th)^2)
            dH_diff = 0.5*(W1_tip - 1.75*W_th)^2;
        end
        % choke losses
        M1_mean_rel = W1_mean/sqrt(gamma*R*T1);
        A_th_star = M1_mean_rel*(A1/N_bl - t)*cos(beta1_geom_mean)/(1+(gamma-1)*M1_mean_rel^2/2)^((gamma+1)/2*(gamma-1))*(1+(gamma-1)/2)^((gamma+1)/2*(gamma-1));
        C_r = sqrt((A1/N_bl-t)*cos(beta1_geom_mean)/A_th);
        if C_r > 1-((A1/N_bl-t)*cos(beta1_geom_mean)/A_th -1)^2
            C_r = 1-((A1/N_bl-t)*cos(beta1_geom_mean)/A_th -1)^2;
        end
        X = 11 - 10*(C_r*A_th)/A_th_star;
        if X<=0
            dH_choke = 0;
        else
            dH_choke = 1/2*W1_mean^2*(0.05*X+X^7);
        end
        %%% PARASSITIC
        % disk friction  (Daily and Nece)
        visc_din2 = mu_NH3(T2); %sutherland
        Re_df = U2*R2/visc_din2*rho2;
        if Re_df < 3e5
            f_df = 2.67/Re_df^0.5;
        else
            f_df = 0.0622/Re_df^0.2;
        end
        dH_disk = f_df * (rho1 + rho2)*R2^2*U2^3/(8*m_dot);
        % recirculation (Coppage)
        D_factor = 1 - W2/W1_tip + 0.75*L_eul*W2/((N_bl/pi*(1-D1_t_opt/D2)+2*D1_t_opt/D2)*W1_tip*U2^2);
        dH_rec = 0.02*sqrt(tan(alpha_2))*D_factor^2*U2^2;
        % leakage (Jansen)
        dH_leak = 0.6 * eps/b2*V2*sqrt(4*pi/(b2*N_bl)*((R1_t_opt-R1_h)/(R2-R1_t_opt))/(1+rho2/rho1)*V2_tg*V1);
        dH_tot_internal =  dH_inc + dH_fr + dH_bl + dH_mix + dH_cl + dH_diff + dH_choke;
        dH_tot_parassitic = dH_disk + dH_rec + dH_leak;
        dH_tot = dH_tot_internal + dH_tot_parassitic;
        eta_new = (L_req - dH_tot) / L_req;
        err_loc = (eta_new - eta_tt_guess)/eta_tt_guess;
        eta_tt_guess = eta_tt_guess + 0.5 * (eta_new - eta_tt_guess);
    end
    eta_sens(i) = eta_tt_guess;
    b2_sens(i) = b2_curr;
end

%% PLOT
figure
plot(D2_vec*1000, eta_sens, 'b-')
grid on
hold on
xline(D2*1000, '--k', 'Balje Reference')
[max_eta, idx] = max(eta_sens);
plot(D2_vec(idx)*1000, max_eta, 'ro', 'MarkerFaceColor','r')
xlabel('D2 [mm]')
ylabel('\eta_{tt}')
title('D2 sensitivity analysis')
legend('Sensitivity','Reference','Optimum','Location','Best')

figure
plot(D2_vec*1000, b2_sens*1000, 'r-', 'LineWidth', 2)
grid on
hold on
plot(D2_vec(idx)*1000, b2_sens(idx)*1000, 'ko', 'MarkerFaceColor','k')
xlabel('D2 [mm]')
ylabel('b2 [mm]')
title('Blade Height Constraint');