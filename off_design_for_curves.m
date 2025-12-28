%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     OFF-DESIGN PERFORMANCE ANALYSIS                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prendo design point
geom.D2 = D2; geom.R2 = D2/2; geom.b2 = b2;
geom.N_bl = N_bl; geom.t = t;
geom.D1_m = D1_m; geom.b1 = b1; 
geom.R1_m = R1_m; geom.R1_t = R1_t_opt; geom.R1_h = R1_h;
geom.beta1_metal = atan(1 - t*N_bl/(2*pi*R1_m)*tan(beta1_mean)); 
geom.beta2_metal = beta2_geom; 
geom.D3 = D3; geom.b3 = b3; geom.R3 = D3/2;
geom.alpha2_set = alpha_2; % Design setting

% costruisco vettore portate
m_dot_design = m_dot;
n_p = 50;
m_dot_vec = linspace(m_dot_design*0.5, m_dot_design*1.4, n_p);


beta_vec = zeros(1, n_p);
eta_vec  = zeros(1, n_p);
surge_vec = zeros(1, n_p);
choke_vec = zeros(1, n_p);
Omega = omega;

% for per grafico
for i = 1:length(m_dot_vec)
    m_curr = m_dot_vec(i);
    

    A1_curr = pi * geom.D1_m * geom.b1;
    V1_curr = m_curr / (rho1 * A1_curr); 
    U1_mean = Omega * geom.R1_m;
    W1_mean_curr = sqrt(U1_mean^2 + V1_curr^2);
    beta1_flow = atan(-U1_mean / V1_curr);

    if abs(geom.beta1_metal - beta1_flow)*180/pi > 12
        surge_vec(i) = 1;
    end
    
    % ora come main
    rho2_loop = rho1 * 2; 
    eta_tt_loop = 0.85; 
    err_eta = 1;
    iter_imp = 0;
    
    while abs(err_eta) > 1e-5
        iter_imp = iter_imp + 1;
        
        V2_m_curr = m_curr / (rho2_loop * pi * geom.D2 * geom.b2);
        U2_curr = Omega * geom.R2;
        mu_slip = 1 - 0.63*pi/geom.N_bl;
        
        % Off-Design: V_tg determinato da beta_metal
        V2_tg_curr = mu_slip * U2_curr + V2_m_curr * tan(geom.beta2_metal);
        
        V2_curr = sqrt(V2_m_curr^2 + V2_tg_curr^2);
        W2_tg_curr = V2_tg_curr - U2_curr;
        W2_curr = sqrt(V2_m_curr^2 + W2_tg_curr^2);
        beta2_flow = atan(W2_tg_curr/V2_m_curr);
        
    
        L_curr = U2_curr * V2_tg_curr; 
        Tt2_curr = Tt1 + L_curr/cp;
        T2_curr = Tt2_curr - V2_curr^2/(2*cp);
        
      
        D2=geom.D2; R2=geom.R2; b2=geom.b2; N_bl=geom.N_bl; t=geom.t;
        U2=U2_curr; V2_tg=V2_tg_curr; V2_meridional=V2_m_curr; V2=V2_curr;
        W2=W2_curr; beta2=beta2_flow; T2=T2_curr; L=L_curr;
        m_dot=m_curr; rho2=rho2_loop;
        
        % Inlet Corrente
        V1=V1_curr; W1_mean=W1_mean_curr; 
        % Approssimazione geometrica per perdite
        D1_t_opt=geom.R1_t*2; D1_h=geom.R1_h*2; D1_m=geom.D1_m; b1=geom.b1;
        R1_t_opt=geom.R1_t; R1_h=geom.R1_h; R1_m=geom.R1_m;
        
        % Angoli e Parametri
        beta1_mean = beta1_flow; % Flow angle corrente
        beta1_tip = atan((-Omega*geom.R1_t)/V1);
        beta1_hub = atan((-Omega*geom.R1_h)/V1);
        
        % Angoli Metal (Fissi) per Incidence/Diffusion
        beta1_geom_mean = geom.beta1_metal;
        beta1_geom_tip  = atan(1-t*N_bl/(2*pi*R1_t_opt)*tan(beta1_tip));
        
        Q_in_curr = m_curr / rho1;
        phi = 4 * Q_in_curr/(pi*U2*D2^2);
        L_eul = L;
        
        dtheta = log(D2/D1_m) * tan(beta2);
        c = R2 * dtheta/(sin(beta2));
        s = c/(1/0.4); 
        
        W1_tip_t = -Omega*R1_t_opt; 
        W1_tip = sqrt(W1_tip_t^2 + V1^2);
        W1_mean_t = -Omega*R1_m;
        
        
        dH_inc = 0.4*(W1_mean - V1/cos(beta1_mean))^2;
       
        visc_din1 = mu_NH3(T1); 
        D_h = D2*cos(beta2) / (( N_bl/pi + D2*cos(beta2)/b2 ) + ( 0.5*(D1_t_opt/D2 + D1_h/D2) * ((cos(beta1_tip) + cos(beta1_hub))/2) )/( N_bl/pi + ((D1_t_opt + D1_h)/(D1_t_opt - D1_h)) * ((cos(beta1_tip) + cos(beta1_hub))/2) ) );
        Re_f = U2 * D_h/visc_din1*rho1;
        Cf = 0.0412 * Re_f ^(-0.1925);
        Lz = 0.08 + 3.16*phi; 
        lb = pi/8 *(D2- 1/2 * (D1_t_opt+D1_h) - b2 + 2*Lz)*(2/((cos(beta1_tip)+cos(beta1_hub))/2 + (cos(beta2))));
        W12_av = 1/4 * (W1_mean + W1_tip + 2*W2); % Approx mean/tip mix
        dH_fr = (2 * Cf * W12_av^2 * lb)/(D_h);
        
        eps = 5e-4;
        dH_cl = 0.6 * eps / b2 * V2_tg * ( 4*pi / (b2 * N_bl) * (R1_t_opt^2 - R1_h^2) / (R2 - R1_t_opt) *V2_tg * V1 / (1 + rho2 / rho1) )^(1/2);
        
        dW = 2*pi*D2*V2_tg/(N_bl*lb);
        dH_bl = 1/48*dW^2;
        
        W_max = 0.5*(W1_mean + W2 + dW);
        D_eq = W_max/W2;
        if D_eq <= 2
            W_sep = W2;
        else
            W_sep = W2*D_eq*0.5; 
        end
        t_te = 2e-3; 
        eps2 = 1- (N_bl * t_te)/(pi*D2);
        A2 = pi*b2*D2*eps2;
        W_out = sqrt((V2_meridional*A2/(pi*D2*b2))^2+W1_mean_t^2);
        dH_mix = 1/2 * (W_sep - W_out)^2;
        
        A_th_geom = s * cos(beta1_geom_mean)*b1;
        A_th = 0.97 * A_th_geom;
        W_th = m_dot/N_bl/rho1/A_th;
        dH_diff = 0.4*(W1_mean - W_th)^2;
        if (W1_tip/W_th > 1.75) && (dH_diff < 0.5*(W1_tip - 1.75*W_th)^2)
            dH_diff = 0.5*(W1_tip - 1.75*W_th)^2;
        end
        
        M1_mean_rel = W1_mean/sqrt(gamma*R*T1);
        A_th_star = M1_mean_rel*(A1/N_bl - t)*cos(beta1_geom_mean)/(1+(gamma-1)*M1_mean_rel^2/2)^((gamma+1)/2*(gamma-1))*(1+(gamma-1)/2)^((gamma+1)/2*(gamma-1));
        C_r = sqrt((A1/N_bl-t)*cos(beta1_geom_mean)/A_th);
        if C_r > 1-((A1/N_bl-t)*cos(beta1_geom_mean)/A_th -1)^2
            C_r = 1-((A1/N_bl-t)*cos(beta1_geom_mean)/A_th -1)^2;
        end
        X = 11 - 10*(C_r*A_th)/A_th_star;
        if X<=0, dH_choke = 0; else, dH_choke = 1/2*W1_mean^2*(0.05*X+X^7); end
        
        visc_din2 = mu_NH3(T2); 
        Re_df = U2*R2/visc_din2*rho2;
        if Re_df < 3e5, f_df = 2.67/Re_df^0.5; else, f_df = 0.0622/Re_df^0.2; end
        dH_disk = f_df * (rho1 + rho2)*R2^2*U2^3/(8*m_dot);
        
        D_factor = 1 - W2/W1_tip + 0.75*L_eul*W2/((N_bl/pi*(1-D1_t_opt/D2)+2*D1_t_opt/D2)*W1_tip*U2^2);
        dH_rec = 0.02*sqrt(tan(geom.alpha2_set))*D_factor^2*U2^2;
        dH_leak = 0.6 * eps/b2*V2*sqrt(4*pi/(b2*N_bl)*((R1_t_opt-R1_h)/(R2-R1_t_opt))/(1+rho2/rho1)*V2_tg*V1);
        
        dH_tot = dH_inc + dH_fr + dH_bl + dH_mix + dH_cl + dH_diff + dH_choke + dH_disk + dH_rec + dH_leak;
        
        % Check Choke Flag
        if dH_choke > 100 || M1_mean_rel > 1.0
            choke_vec(i) = 1;
        end
        
        % Update
        eta_new = (L_curr - dH_tot) / L_curr;
        Pt2_curr = Pt1 * (1 + eta_new * L_curr/(cp*Tt1))^(gamma/(gamma-1));
        M2_curr = V2_curr / sqrt(gamma*R*T2_curr);
        P2_curr = Pt2_curr / (1 + (gamma-1)/2 * M2_curr^2)^(gamma/(gamma-1));
        rho2_new = P2_curr / (R * T2_curr);
        
        err_eta = (eta_new - eta_tt_loop)/eta_tt_loop;
        eta_tt_loop = eta_tt_loop + 0.5 * (eta_new - eta_tt_loop);
        rho2_loop = 0.5 * rho2_loop + 0.5 * rho2_new;
    end
    
    % --- VANELESS DIFFUSER LOOP (Identico al Main) ---
    D3 = geom.D3; b3 = geom.b3; R3 = geom.R3;
    Tt3 = Tt2_curr; Pt2 = Pt2_curr;
    rho3 = rho2_new; V3 = V2_curr; 
    errV3 = 1; errrho3 = 1;
    k_calib = 0.01;
    iter_vd = 0;
    
    while (abs(errrho3) > 1e-5 || abs(errV3) > 1e-5) && iter_vd < 50
        iter_vd = iter_vd + 1;
        rho_avg = 0.5*(rho2_new + rho3);
        V_avg = 0.5 * (V2_curr + V3);
        D_avg = 0.5 * (geom.D2 + D3);
        
        T_avg_diff = 0.5 * (T2_curr + (Tt3 - V3^2/(2*cp)));
        visc_din3 = mu_NH3(T_avg_diff);
        
        Re_avg = rho_avg * V_avg * D_avg/visc_din3;
        Cf_vaned = k_calib*(1.8e5/Re_avg)^0.2;
        
        V3_tg = V2_tg_curr/(D3/geom.D2 + 0.5*pi*Cf_vaned*rho2_new*V2_tg_curr*D3*(D3-geom.D2)/m_curr);
        V3_meridional = (V2_m_curr*rho2_new*pi*geom.D2*geom.b2)/(rho3*pi*D3*b3);
        
        V3_calc = sqrt(V3_meridional^2 + V3_tg^2);
        errV3 = (V3_calc - V3)/V3;
        V3 = V3 + 0.6 * (V3_calc - V3);
        
        dH_t_VD = Cf_vaned * V1^2 * R2 *(1-(R2/R3)^1.5)/(1.5*b2*cos(geom.alpha2_set));
        Tt3_is = Tt3 - dH_t_VD/cp;
        Pt3 = Pt2 * (Tt3_is/Tt2_curr)^(gamma/(gamma-1));
        
        T3 = Tt3 - 1/(2*cp)*V3^2;
        P3 = Pt3/(1+(gamma-1)/2*(V3/sqrt(gamma*R*T3))^2)^(gamma/(gamma-1));
        
        rho3_calc = P3/(T3*R);
        errrho3 = (rho3_calc - rho3)/rho3;
        rho3 = rho3 + 0.6 * (rho3_calc - rho3);
    end
    
    % --- VANED DIFFUSER (Performance Off-Design) ---
    % Nota: Usiamo C_p o eta_d fissi dal main per ottenere Pt4
    C_p = 0.8; 
    P4 = P3 + C_p*(Pt3 - P3);
    T4_is = T3*(P4/P3)^((gamma-1)/gamma);
    M3 = V3/sqrt(gamma*R*T3); % Mach ingresso vaned
    
    % Pressione Totale Uscita Stage
    % Approx: Pt4 dipende da P4 e Mach uscita (che dipende da S4 fissa)
    % Usiamo relazione isentropica su P4
    Pt4_curr = P4*(1 + (gamma-1)/2 * (M3*0.6)^2)^(gamma/(gamma-1)); % Stima M4 < M3
    
    % Global Stage Efficiency
    L_is_stage = cp*Tt1*((Pt4_curr/Pt1)^((gamma-1)/gamma)-1);
    eta_stage = L_is_stage / L_curr;
    
    beta_vec(i) = Pt4_curr / Pt1;
    eta_vec(i) = eta_stage;
end

% 4. Plot
figure('Color','w')
subplot(2,1,1)
plot(m_dot_vec, beta_vec, 'b-o'); grid on; hold on;
xlabel('Mass Flow [kg/s]'); ylabel('Stage PR \beta_{tt}');
title('Stage Pressure Ratio');
idx_surge = find(surge_vec==1, 1, 'last');
if ~isempty(idx_surge), xline(m_dot_vec(idx_surge), '--r', 'Surge'); end
idx_choke = find(choke_vec==1, 1, 'first');
if ~isempty(idx_choke), xline(m_dot_vec(idx_choke), '--k', 'Choke'); end

subplot(2,1,2)
plot(m_dot_vec, eta_vec, 'g-o'); grid on; hold on;
xlabel('Mass Flow [kg/s]'); ylabel('Stage Efficiency');
title('Stage Efficiency');