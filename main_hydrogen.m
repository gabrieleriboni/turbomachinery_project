%% hydrogen.m
clear
close all
clc

%% Load geometry/design point
S = load('geom_design.mat');
geom   = S.geom;
design = S.design;

%% 1. Gas properties HYDROGEN
% Valori per l'H2
cp_exp_H2 = 14.310; % kJ/(kg.K) (oppure in moli e convertito)
cv_exp_H2 = 10.184; % kJ/(kg.K)
Mm_H2     = 2.016;  % g/mol

cp = cp_exp_H2 * 1e3; % J/(kg.K)
cv = cv_exp_H2 * 1e3; % J/(kg.K)
gamma = cp/cv;
R = cp - cv;          % circa 4126 J/(kg.K)

%% 2. Sutherland's Law for Hydrogen
mu_0_H2 = 8.76e-06;
T_0_H2  = 293.85;
S_H2    = 72;
mu_gas = @(T) mu_0_H2 * (T./T_0_H2).^(1.5) .* (T_0_H2 + S_H2) ./ (T + S_H2);

%% 3. Operating condition (fixed speed + inlet totals)
Pt1   = geom.Pt1;
Tt1   = geom.Tt1;
omega = geom.omega;

%% 4. Mass-flow sweep (CORRETTO PER SIMILITUDINE H2)
% Per l'ammoniaca avevi R_NH3 circa 528.4 J/kgK (da cp_exp=0.037, cv_exp=0.028, Mm=17.03)
R_NH3 = ((0.037 - 0.028)/17.03) * 1e6; 

% La portata di design volumetrica deve rimanere costante.
% Quindi m_dot_H2 = m_dot_NH3 * (rho_H2 / rho_NH3) = m_dot_NH3 * (R_NH3 / R_H2)
m_dot_des_H2 = geom.m_dot_design * (R_NH3 / R); 

fprintf('\n>>> PORTATA EQUIVALENTE CALCOLATA PER IDROGENO: %.4f kg/s <<<\n\n', m_dot_des_H2);

% Sweep centrato sulla nuova portata dell'H2 (circa 0.38 kg/s)
m_dot_vec = sort(unique([linspace(0.6*m_dot_des_H2, 1.4*m_dot_des_H2, 81), m_dot_des_H2]));

% Preallocate
beta_tt_vec   = nan(size(m_dot_vec));
eta_c_vec     = nan(size(m_dot_vec));
eta_stage_vec = nan(size(m_dot_vec));
X_imp_vec     = nan(size(m_dot_vec));
X_vd_vec      = nan(size(m_dot_vec));
choked_vec    = false(size(m_dot_vec));
converged_vec = false(size(m_dot_vec));

fprintf('OFF-DESIGN sweep started...\n');
for ii = 1:numel(m_dot_vec)
    m_dot = m_dot_vec(ii);
    res = one_point_fixed_geom(m_dot, geom, Pt1, Tt1, omega, cp, gamma, R, mu_gas);
    
    beta_tt_vec(ii)   = res.beta_tt;
    eta_c_vec(ii)     = res.eta_c;
    eta_stage_vec(ii) = res.eta_stage_tt;
    X_imp_vec(ii)     = res.X_imp;
    X_vd_vec(ii)      = res.X_vd;
    choked_vec(ii)    = res.is_choked;
    converged_vec(ii) = res.is_converged;
    
    fprintf('m_dot=%.4f | beta_tt=%.4f | eta=%.4f | Ximp=%.3f | Xvd=%.3f | Mmax=%.3f | choked=%d | conv=%d\n', ...
        m_dot, res.beta_tt, res.eta_c, res.X_imp, res.X_vd, res.M_max, res.is_choked, res.is_converged);
end

%% Keep only good points
ok = converged_vec ...
     & isfinite(beta_tt_vec) & isreal(beta_tt_vec) & beta_tt_vec > 1 & beta_tt_vec < 10 ...
     & isfinite(eta_c_vec)   & isreal(eta_c_vec)   & eta_c_vec > 0   & eta_c_vec < 1.2;

%% Dimensionless Mass Flow Parameter
rho_t1 = Pt1 / (R * Tt1);
a_t1 = sqrt(gamma * R * Tt1);
m_nd_vec = m_dot_vec ./ (rho_t1 * a_t1 * geom.D2^2);

m_nd_ok   = m_nd_vec(ok);
b_ok   = beta_tt_vec(ok);
eta_ok = eta_c_vec(ok);

%% Diagnostics
% We mark choke when any station Mach number (M1, M2, M3, M4, M6) >= 1.
choke_idx = find(choked_vec, 1, 'first');
if ~isempty(choke_idx)
    m_nd_choke = m_nd_vec(choke_idx);
else
    m_nd_choke = nan;
end

[~, i_peak] = max(b_ok);
if ~isempty(i_peak)
    m_nd_surge = m_nd_ok(i_peak);
else
    m_nd_surge = nan;
end

%% Calcolo del vero punto di design per H2 (Picco di efficienza)
if ~isempty(eta_ok)
    [max_eta_H2, idx_des_H2] = max(eta_ok);
    m_nd_opt_H2 = m_nd_ok(idx_des_H2);
    m_dot_opt_H2 = m_nd_opt_H2 * (rho_t1 * a_t1 * geom.D2^2);
    beta_opt_H2 = b_ok(idx_des_H2);
    
    [~, idx_sim] = min(abs(m_dot_vec - m_dot_des_H2));
    m_nd_sim_H2 = m_nd_vec(idx_sim);
    beta_sim_H2 = beta_tt_vec(idx_sim);
    eta_sim_H2 = eta_c_vec(idx_sim);
    
    fprintf('\n==========================================\n');
    fprintf('RISULTATI VERIFICA IDROGENO (A PARI RPM):\n');
    fprintf('Portata ottimale H2 (Picco): %.4f kg/s (m_nd = %.4f)\n', m_dot_opt_H2, m_nd_opt_H2);
    fprintf('Efficienza massima H2:       %.4f\n', max_eta_H2);
    fprintf('Rapporto Compressione Picco: %.4f\n', beta_opt_H2);
    fprintf('------------------------------------------\n');
    fprintf('PUNTO DI SIMILITUDINE GEOMETRICA:\n');
    fprintf('Portata similitudine:        %.4f kg/s (m_nd = %.4f)\n', m_dot_des_H2, m_nd_sim_H2);
    fprintf('Efficienza similitudine:     %.4f\n', eta_sim_H2);
    fprintf('Rapporto Compressione sim:   %.4f\n', beta_sim_H2);
    fprintf('==========================================\n\n');
end

%% Plots
% Plot: flow - pressure ratio
figure
plot(m_nd_vec, beta_tt_vec, '.-')
grid on
xlabel('Dimensionless Mass Flow $\dot{m}_{nd} = \frac{\dot{m}}{\rho_{t1} a_{t1} D_2^2}$ [-]', 'Interpreter', 'latex')
ylabel('$\beta_{tt} = P_{t,out}/P_{t,in}$ [-]', 'Interpreter', 'latex')
title('\textbf{H2 Compressor Curve: \boldmath$\dot{m}_{nd} - \beta_{tt}$}', 'Interpreter', 'latex')
hold on

if exist('m_nd_opt_H2', 'var')
    plot(m_nd_opt_H2, beta_opt_H2, 'o', 'Color', [0.7, 0, 0], 'MarkerSize', 8, 'LineWidth', 1.5)
    text(m_nd_opt_H2, beta_opt_H2, sprintf('\\quad Peak ($\\dot{m}_{nd}=%.4f$, %.3f)', m_nd_opt_H2, beta_opt_H2), 'Interpreter', 'latex', 'Color', [0.7, 0, 0], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
    
    plot(m_nd_sim_H2, beta_sim_H2, 's', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'LineWidth', 1.5)
    text(m_nd_sim_H2, beta_sim_H2, sprintf('\\quad Similitude ($\\dot{m}_{nd}=%.4f$, %.3f)', m_nd_sim_H2, beta_sim_H2), 'Interpreter', 'latex')
end

if ~isnan(m_nd_surge), xline(m_nd_surge, '--', 'Surge (peak \beta_{tt})'); end
if ~isnan(m_nd_choke), xline(m_nd_choke, '--', 'Choke (M \geq 1)'); end

% Plot: flow - efficiency
figure
plot(m_nd_vec, eta_c_vec, '.-')
grid on
xlabel('Dimensionless Mass Flow $\dot{m}_{nd} = \frac{\dot{m}}{\rho_{t1} a_{t1} D_2^2}$ [-]', 'Interpreter', 'latex')
ylabel('$\eta$ [-]', 'Interpreter', 'latex')
title('\textbf{H2 Compressor Curve: \boldmath$\dot{m}_{nd} - \eta$}', 'Interpreter', 'latex')
hold on

if exist('m_nd_opt_H2', 'var')
    plot(m_nd_opt_H2, max_eta_H2, 'o', 'Color', [0.7, 0, 0], 'MarkerSize', 8, 'LineWidth', 1.5)
    text(m_nd_opt_H2, max_eta_H2, sprintf('\\quad Peak ($\\dot{m}_{nd}=%.4f$, %.4f)', m_nd_opt_H2, max_eta_H2), 'Interpreter', 'latex', 'Color', [0.7, 0, 0], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
    
    plot(m_nd_sim_H2, eta_sim_H2, 's', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'LineWidth', 1.5)
    text(m_nd_sim_H2, eta_sim_H2, sprintf('\\quad Similitude ($\\dot{m}_{nd}=%.4f$, %.4f)', m_nd_sim_H2, eta_sim_H2), 'Interpreter', 'latex')
end

if ~isnan(m_nd_surge), xline(m_nd_surge, '--', 'Surge (peak \beta_{tt})'); end
if ~isnan(m_nd_choke), xline(m_nd_choke, '--', 'Choke (M \geq 1)'); end

%% ====================== LOCAL FUNCTION ======================
function res = one_point_fixed_geom(m_dot, geom, Pt1, Tt1, omega, cp, gamma, R, mu_gas)
    res = struct();
    res.is_converged = true;
    res.is_choked    = false;
    res.beta_tt      = NaN;
    res.eta_c        = NaN;
    res.eta_stage_tt = NaN;
    res.X_imp        = NaN;
    res.X_vd         = NaN;
    res.M_max        = NaN;
    
    D1_h = geom.D1_h; R1_h = geom.R1_h;
    R1_t_opt = geom.R1_t; D1_t_opt = geom.D1_t;
    R1_m = geom.R1_m; D1_m = geom.D1_m;
    b1 = geom.b1; A1 = geom.A1;
    D2 = geom.D2; R2 = geom.R2; b2 = geom.b2; U2 = omega * R2;
    N_bl = geom.N_bl; t = geom.t; t_te = geom.t_te; eps = geom.eps;
    beta2_geom_fix = geom.beta2_geom;
    R3 = geom.R3; D3 = geom.D3; b3_fix = geom.b3;
    N_bl_VD = geom.N_bl_VD;
    R4 = geom.R4; D4 = geom.D4; b4 = geom.b4;
    alpha3_geom = geom.alpha3_geom; alpha4_geom = geom.alpha4_geom;
    W3 = geom.W3; W4 = geom.W4; Lb = geom.Lb; theta_c = geom.theta_c;
    R5 = geom.R5; D5 = geom.D5; b5 = geom.b5;
    SP = geom.SP; K_friction = geom.K_friction;
    r6 = geom.r6; A6 = geom.A6; D_throat = geom.D_throat; AR_cone = geom.AR_cone;
    
    rho_t1 = Pt1/(R*Tt1);
    
    %% 1) Impeller inlet
    rho1 = rho_t1;
    err_new = 1; relaxation = 0.8; it = 0;
    while abs(err_new) > 1e-5 && it < 200
        rho = rho1;
        cm1 = m_dot ./ (rho * pi * (R1_t_opt^2 - R1_h^2));
        T1 = Tt1 - cm1^2/(2*cp);
        if ~isreal(T1) || T1 <= 1
            res.is_choked = true; res.is_converged = false; return
        end
        P1 = Pt1 * (T1/Tt1)^(gamma/(gamma-1));
        rho_new = P1/(R*T1);
        rho_new = rho + relaxation * (rho_new - rho);
        err_new = abs(rho_new - rho)/rho;
        rho1 = rho_new;
        it = it + 1;
    end
    if it >= 200, res.is_converged = false; return; end
    
    U1_tip = omega * R1_t_opt; U1_mean = omega * R1_m; U1_hub = omega * R1_h;
    W1_tip_t = -U1_tip; W1_mean_t = -U1_mean; W1_hub_t = -U1_hub;
    W1_meridional = m_dot/rho1/pi/(R1_t_opt^2-R1_h^2);
    W1_tip = sqrt(W1_tip_t^2 + W1_meridional^2);
    W1_mean = sqrt(W1_mean_t^2 + W1_meridional^2);
    W1_hub = sqrt(W1_hub_t^2 + W1_meridional^2);
    beta1_tip = atan(W1_tip_t/W1_meridional);
    beta1_mean = atan(W1_mean_t/W1_meridional);
    beta1_hub = atan(W1_hub_t/W1_meridional);
    V1 = W1_meridional;
    M1 = V1/sqrt(gamma*R*T1);
    
    %% 2) Impeller outlet
    eta_tt = 0.90; rho2 = rho1;
    err_eta_new = 1; err_rho2_new = 1; relaxation2 = 0.8; relaxation_rho2 = 0.8; k = 0;
    X_imp = -inf;
    while (abs(err_eta_new) > 1e-5 || abs(err_rho2_new) > 1e-5) && k < 300
        V2_meridional = m_dot/(rho2*pi*D2*b2);
        mu = 1 - 0.63*pi/N_bl;
        V2_tg = mu*U2 + V2_meridional*tan(beta2_geom_fix);
        V2 = sqrt(V2_meridional^2 + V2_tg^2);
        alpha2 = atan(V2_tg/V2_meridional);
        W2_tg = V2_tg - U2; W2_meridional = V2_meridional;
        W2 = sqrt(W2_meridional^2 + W2_tg^2);
        beta2 = atan(W2_tg/V2_meridional);
        
        L_eul = U2*V2_tg; L = L_eul;
        Tt2 = L/cp + Tt1; T2 = Tt2 - V2^2/(2*cp);
        if ~isreal(T2) || T2 <= 1
            res.is_choked = true; res.is_converged = false; return
        end
        M2 = V2/sqrt(gamma*R*T2);
        
        basePt2 = 1 + (eta_tt*L)/(cp*Tt1);
        if ~isreal(basePt2) || basePt2 <= 0
            res.is_converged = false; return
        end
        Pt2 = Pt1 * (basePt2)^(gamma/(gamma - 1));
        P2 = Pt2/(1 + (gamma-1)/2 * M2^2)^(gamma/(gamma - 1));
        rho2_new = P2/(R*T2);
        rho2_new_rel = rho2 + relaxation_rho2*(rho2_new - rho2);
        err_rho2_new = abs(rho2_new_rel - rho2)/rho2;
        rho2 = rho2_new_rel;
        
        dH_inc = 0.4*(W1_mean - V1/cos(geom.beta1_geom_mean))^2;
        visc_din1 = mu_gas(Tt1 - V1^2/(2*cp));
        term_den1 = N_bl/pi + D2*cos(beta2)/b2;
        term_num = 0.5*(D1_t_opt/D2 + D1_h/D2) * ((cos(beta1_tip) + cos(beta1_hub))/2);
        term_den2 = N_bl/pi + ((D1_t_opt + D1_h)/(D1_t_opt - D1_h)) * ((cos(beta1_tip) + cos(beta1_hub))/2);
        D_h = D2*(cos(beta2) / term_den1 + term_num/term_den2);
        Re_f = U2 * D_h/visc_din1*rho_t1;
        Cf = 0.0412 * Re_f ^(-0.1925);
        phi = 4 * (m_dot/rho_t1)/(pi*U2*D2^2);
        Lz = 0.08 + 3.16*phi;
        lb = pi/8 *(D2- 1/2 * (D1_t_opt+D1_h) - b2 + 2*Lz)*(2/((cos(beta1_tip)+cos(beta1_hub))/2 + cos(beta2)));
        W12_av = 1/4 * (W1_hub + W1_tip + 2*W2);
        dH_fr = (2 * Cf * W12_av^2 * lb)/D_h;
        dH_cl = 0.6 * eps / b2 * V2_tg * sqrt(4*pi / (b2 * N_bl) * (R1_t_opt^2 - R1_h^2) / (R2 - R1_t_opt) *V2_tg * V1 / (1 + rho2 / rho1));
        dW = 2*pi*D2*V2_tg/(N_bl*lb); dH_bl = 1/48*dW^2;
        W_max = 0.5*(W1_mean + W2 + dW); D_eq = W_max/W2;
        if D_eq <= 2, W_sep = W2; else, W_sep = W2*D_eq*0.5; end
        eps2 = 1- (N_bl * t_te)/(pi*D2);
        A2_annulus = pi*b2*D2; A2_passages = A2_annulus * eps2;
        W_out = sqrt(V2_meridional^2 + W2_tg^2); % V2_meridional is already mixed-out
        dH_mix = 1/2 * (W_sep - W_out)^2;
        A_th = 0.97 * geom.A_th_imp_geom;
        W_th = m_dot/N_bl/rho1/A_th; dH_diff = 0.4*(W1_mean - W_th)^2;
        if (W1_tip/W_th > 1.75) && (dH_diff < 0.5*(W1_tip - 1.75*W_th)^2)
            dH_diff = 0.5*(W1_tip - 1.75*W_th)^2;
        end
        M1_mean_rel = W1_mean/sqrt(gamma*R*(Tt1 - V1^2/(2*cp)));
        A_th_star = M1_mean_rel*(A1/N_bl - t)*cos(geom.beta1_geom_mean) / ( (2/(gamma+1))*(1+(gamma-1)/2*M1_mean_rel^2) )^((gamma+1)/(2*(gamma-1)));
        C_r = sqrt((A1/N_bl-t)*cos(geom.beta1_geom_mean)/A_th);
        if C_r > 1-((A1/N_bl-t)*cos(geom.beta1_geom_mean)/A_th -1)^2
            C_r = 1-((A1/N_bl-t)*cos(geom.beta1_geom_mean)/A_th -1)^2;
        end
        X_imp = 11 - 10*(C_r*A_th)/A_th_star;
        if X_imp<=0, dH_choke = 0; else, dH_choke = 1/2*W1_mean^2*(0.05*X_imp+X_imp^7); end
        visc_din2 = mu_gas(T2);
        Re_df = U2*R2/visc_din2*rho2;
        if Re_df < 3e5, f_df = 2.67/Re_df^0.5; else, f_df = 0.0622/Re_df^0.2; end
        dH_disk = f_df * (rho1 + rho2)*R2^2*U2^3/(8*m_dot);
        D_factor = 1 - W2/W1_tip + 0.75*L_eul*W2/((N_bl/pi*(1-D1_t_opt/D2)+2*D1_t_opt/D2)*W1_tip*U2^2);
        dH_rec = 0.02*sqrt(tan(alpha2))*D_factor^2*U2^2;
        dH_leak = 0.6 * eps/b2*V2*sqrt(4*pi/(b2*N_bl)*((R1_t_opt-R1_h)/(R2-R1_t_opt))/(1+rho2/rho1)*V2_tg*V1);
        
        dH_tot_internal =  dH_inc + dH_fr + dH_bl + dH_mix + dH_cl + dH_diff + dH_choke;
        dH_tot_parassitic = dH_disk + dH_rec + dH_leak;
        eta_tt_new = (L_eul - dH_tot_internal)/(L_eul + dH_tot_parassitic);
        eta_tt_new = eta_tt + relaxation2 * (eta_tt_new - eta_tt);
        err_eta_new = abs(eta_tt_new - eta_tt)/eta_tt;
        eta_tt = eta_tt_new;
        k = k + 1;
    end
    if k >= 300, res.is_converged = false; return; end
    
    %% 3) Vaneless diffuser
    b3 = b3_fix; Tt3 = Tt2; visc_din3 = visc_din2; rho3 = rho2; V3 = V2;
    k_vld = 0.01; errV3 = 1; errrho3 = 1; j = 0; relaxation3 = 0.8;
    while (abs(errrho3)>1e-5 || abs(errV3)>1e-5) && j < 200
        rho_avg = 0.5*(rho2 + rho3); V_avg = 0.5*(V2 + V3); D_avg = 0.5*(D2+D3);
        Re_avg = rho_avg * V_avg * D_avg/visc_din3;
        Cf_vaneless = k_vld*(1.8e5/Re_avg)^0.2;
        alpha3 = atan(b3*tan(alpha2)*rho3/b2/rho2);
        V3_meridional = (V2_meridional*rho2*pi*D2*b2)/(rho3*pi*D3*b3);
        V3_tg = V3_meridional*tan(alpha3); V3_new = V3_meridional/cos(alpha3);
        errV3 = abs(V3_new - V3)/max(V3,1e-9); V3 = V3 + 0.8*(V3_new - V3);
        dH_t_VLD = Cf_vaneless * V2^2 * R2 *(1-(R2/R3)^1.5)/(1.5*b2*cos(alpha2));
        T3 = Tt3-1/(2*cp)*V3^2;
        if ~isreal(T3) || T3 <= 1, res.is_choked = true; res.is_converged = false; return; end    
        M3 = V3/sqrt(gamma*R*T3);
        Tt3_is = Tt3-dH_t_VLD/cp;
        Pt3 = Pt2*(Tt3_is/Tt2)^(gamma/(gamma-1));
        P3 = Pt3/(1+(gamma-1)/2*M3^2)^(gamma/(gamma-1));
        rho3_new = P3/(T3*R);
        errrho3 = abs(rho3_new - rho3)/rho3;
        rho3 = rho3 + relaxation3*(rho3_new - rho3);
        j = j+1;
    end
    if j >= 200, res.is_converged=false; return; end
    
    %% 4) Vaned diffuser
    delta_s = (geom.theta * (0.92*(geom.camber)^2+0.02*alpha4_geom))/(sqrt(geom.sigma)-0.02*geom.theta);
    d_delta_i = exp(((1.5-(90-alpha3_geom)/60)^2-3.3)*geom.sigma);
    alpha4 = alpha4_geom - delta_s - d_delta_i * (alpha3_geom - alpha3);
    rho4 = rho3; err_rho4 = 1; it_vd = 0; X_vd = -inf;
    while err_rho4 > 1e-5 && it_vd < 80
        rho4_old = rho4;
        V4_meridional = m_dot / (rho4_old * N_bl_VD * W4 * b4);
        V4 = V4_meridional / cos(alpha4); V4_tg = V4 * sin(alpha4);
        A_inlet_VD = W3 * b3_fix; cos_alpha_th = cos(alpha3)^2/cos(alpha3_geom);
        A_inlet_geom = 2 * pi * R3 * b3_fix; A_throat_geom = N_bl_VD * W3 * b3_fix;
        cos_alpha_th_geom = A_throat_geom / A_inlet_geom;
        cos_alpha_opt = sqrt(cos(alpha3_geom) * cos_alpha_th_geom);
        V3_opt = V3_meridional / cos_alpha_opt;
        dH_inc_VD = 0.4 * (V3 - V3_opt)^2;
        Dh_in = 2 * (W3 * b3_fix) / (W3 + b3_fix); Dh_out = 2 * (W4 * b4) / (W4 + b4);
        Dh_avg_VD = 0.5 * (Dh_in + Dh_out);
        
        rho_avg_VD = 0.5 * (rho3 + rho4);
        V_avg_VD = 0.5 * (V3 + V4);
        Re_VD = rho_avg_VD * V_avg_VD * Dh_avg_VD / visc_din3;
        Cf_vaned = k_vld*(1.8e5/Re_VD)^0.2;
        
        dH_fr_VD = Cf_vaned * Lb/Dh_avg_VD*((V3+V4)/2)^2;
        Cr_VD = sqrt( (A_inlet_geom * cos(alpha3_geom)) / A_throat_geom );
        M3_abs = V3/sqrt(gamma*R*T3);
        A_th_star = M3_abs * A_throat_geom / ( (2/(gamma+1))*(1+(gamma-1)/2*M3_abs^2) )^((gamma+1)/(2*(gamma-1)));
        X_vd = 11 - 10 * (Cr_VD * A_throat_geom) / A_th_star;
        if X_vd <= 0, dH_choke_VD = 0; else, dH_choke_VD = 1/2 * V3^2 * (0.05 * X_vd + X_vd^7); end
        Cl = 1; C_theta = 1; k1 = 0.2*(1-1/Cl/C_theta); C_theta = max(1, 2*theta_c/0.192);   % theta_c in rad
k2 = (2.88*theta_c)/(pi*C_theta) * (1 - theta_c/(0.192*C_theta));
        Cr = 1/2*V3_meridional*cos(alpha4_geom)/(V4_meridional*cos(alpha3_geom))+1;
        B4 = (k1+k2*(Cr^2-1))*Lb/W4; dH_bl_VD = 1/2*(V4/(1-B4)-V4)^2;
        V_sep = V4/(1+2*C_theta); V4_meridional_mixed = m_dot / (rho4 * 2 * pi * R4 * b4);
        V_out = sqrt(V4_meridional_mixed^2 + V4_tg^2); dH_mix_VD = 0.5 * (V_sep - V_out)^2;
        dH_t_VD = dH_inc_VD + dH_fr_VD + dH_choke_VD + dH_bl_VD + dH_mix_VD;
        eta_int_tt = (L_eul - dH_tot_internal - dH_t_VLD - dH_t_VD)/(L_eul + dH_tot_parassitic);
        Tt4 = Tt2; T4 = Tt4 - V4^2/(2*cp);
        if ~isreal(T4) || T4 <= 1, res.is_choked = true; res.is_converged=false; return; end
        Pt4 = Pt1 * (1 + eta_int_tt * ((Tt2 - Tt1)/Tt1))^(gamma/(gamma-1));
        M4 = V4/sqrt(gamma*R*T4); P4 = Pt4 / (1 + (gamma-1)/2 * M4^2)^(gamma/(gamma-1));
        rho4 = P4/(R*T4); err_rho4 = abs(rho4 - rho4_old)/rho4_old;
        it_vd = it_vd + 1;
    end
    if it_vd >= 80, res.is_converged=false; return; end
    
    %% 5) Second vaneless diffuser
    Tt5 = Tt4; visc_din5 = mu_gas(T4); 
    
    % V4 for the vaneless space based on full annulus mass conservation
    V4_vld_meridional = m_dot / (rho4 * pi * D4 * b4);
    V4_vld_tg = V4_tg; % Tangential velocity is conserved
    V4_vld = sqrt(V4_vld_meridional^2 + V4_vld_tg^2);
    alpha4_vld = atan(V4_vld_tg / V4_vld_meridional);

    rho5 = rho4; V5 = V4_vld;
    kk = 0.01; errV5 = 1; errrho5 = 1; jj = 0; relaxation5 = 0.8;
    while (abs(errrho5)>1e-5 || abs(errV5)>1e-5) && jj < 200
        rho_avg = 0.5*(rho4 + rho5); V_avg = 0.5 * (V4_vld + V5); D_avg = 0.5 * (D4 + D5);
        Re_avg = rho_avg * V_avg * D_avg/visc_din5; Cf2_vaneless = kk*(1.8e5/Re_avg)^0.2;
        alpha5 = atan(rho5/rho4*tan(alpha4_vld));
        V5_meridional = (V4_vld_meridional*rho4*pi*D4*b4)/(rho5*pi*D5*b5);
        V5_tg = V5_meridional*tan(alpha5); V5_new = V5_meridional/cos(alpha5);
        errV5 = abs(V5_new - V5)/max(V5,1e-9); V5 = V5 + 0.8*(V5_new - V5);
        dH_t_2VLD = Cf2_vaneless * V4_vld^2 * R4 *(1-(R4/R5)^1.5)/(1.5*b4*cos(alpha4_vld));
        T5 = Tt5-1/(2*cp)*V5^2;
        if ~isreal(T5) || T5 <= 1, res.is_converged=false; return; end
        M5 = V5/sqrt(gamma*R*T5);
        Tt5_new = Tt5 - dH_t_2VLD/cp; Pt5 = Pt4*(Tt5_new/Tt2)^(gamma/(gamma-1));
        P5 = Pt5/(1+(gamma-1)/2*M5^2)^(gamma/(gamma-1));
        rho5_new = P5/(T5*R); errrho5 = abs(rho5_new - rho5)/rho5;
        rho5 = rho5 + relaxation5*(rho5_new - rho5);
        jj = jj+1;
    end
    if jj >= 200, res.is_converged=false; return; end
    
    %% 6) Volute
    rho6 = rho5; err_vol = 1; iter_vol = 0; relaxation6 = 0.8;
    while err_vol > 1e-5 && iter_vol < 120
        V6 = m_dot/(rho6*A6); Tt6 = Tt5; T6  = Tt6 - V6^2/(2*cp);
        if ~isreal(T6) || T6 <= 1, res.is_converged=false; return; end
        dH_vol_merid = 0.5 * V5_meridional^2;
        if SP >= 1, term_sp = (1 - 1/SP^2); dH_vol_tan = 0.5 * (R5/r6) * V5_tg^2 * term_sp * 0.5;
        else, term_sp = (1 - 1/SP)^2; dH_vol_tan = (R5/r6) * V5_tg^2 * term_sp * 0.5; end
        L_vol = pi * (R5 + r6); d_H = D_throat;
        dH_vol_fric = 2 * K_friction * (L_vol / d_H) * V6^2;
        dH_vol_tot = dH_vol_merid + dH_vol_tan + dH_vol_fric;
        Tt6_is = Tt6 - dH_vol_tot/cp; Pt6 = Pt5 * (Tt6_is/Tt5)^(gamma/(gamma-1));
        M6 = V6/sqrt(gamma*R*T6); P6 = Pt6/(1+(gamma-1)/2*M6^2)^(gamma/(gamma-1));
        rho6_new = P6/(R*T6); err_vol = abs(rho6_new - rho6)/rho6;
        rho6 = rho6 + relaxation6*(rho6_new - rho6);
        iter_vol = iter_vol + 1;
    end
    if iter_vol >= 120, res.is_converged=false; return; end
    
    %% 7) Exit cone
    A7 = A6 * AR_cone; 
    
    rho7 = rho6; err_rho7 = 1; iter_cone = 0;
    while err_rho7 > 1e-5 && iter_cone < 100
        V7 = m_dot/(rho7*A7);
        dH_cone = 0.5 * (V6 - V7)^2; Tt7 = Tt6;
        Tt7_is = Tt7 - dH_cone/cp; Pt7 = Pt6 * (Tt7_is/Tt7)^(gamma/(gamma-1));
        T7 = Tt7 - V7^2/(2*cp);
        if ~isreal(T7) || T7 <= 1, res.is_converged=false; return; end
        M7 = V7/sqrt(gamma*R*T7); P7 = Pt7/(1+(gamma-1)/2*M7^2)^(gamma/(gamma-1));
        rho7_new = P7/(R*T7); err_rho7 = abs(rho7_new - rho7)/rho7;
        rho7 = rho7 + 0.8*(rho7_new - rho7);
        iter_cone = iter_cone + 1;
    end
    if iter_cone >= 100, res.is_converged=false; return; end
    rho7_final = rho7;
    beta_tt = Pt7/Pt1;
    
    L_is = cp*Tt1*(beta_tt^((gamma-1)/gamma)-1);
    eta_c = L_is / L;
    eta_stage_tt = (L_eul - dH_tot_internal - dH_t_VLD - dH_t_VD - dH_t_2VLD - dH_vol_tot - dH_cone)/(L_eul + dH_tot_parassitic);
    
    is_choked = (M1 >= 1) || (M2 >= 1) || (M3 >= 1) || (M4 >= 1) || (M6 >= 1);
    
    res.beta_tt      = beta_tt;
    res.eta_c        = eta_c;
    res.eta_stage_tt = eta_stage_tt;
    res.X_imp        = X_imp;
    res.X_vd         = X_vd;
    res.is_choked    = is_choked;
    res.M_max        = max([M1, M2, M3, M4, M6]);
end