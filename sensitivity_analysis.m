% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                     DETAILED SENSITIVITY ANALYSIS                       %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clc
% 
% n_sens = 50;  
% range_var = 0.20; % +/- 20%
% D2_vec = linspace(D2 * (1 - range_var), D2 * (1 + range_var), n_sens);
% 
% % Vettori output
% eta_sens = zeros(1, n_sens);
% b2_sens  = zeros(1, n_sens);
% beta2_sens = zeros(1, n_sens);
% 
% % Dati fissi dall'Inlet (non cambiano variando solo D2)
% U1_mean_fix = omega * R1_m;
% W1_mean_fix = sqrt(U1_mean_fix^2 + V1^2); % V1 = V_meridional_1
% beta1_mean_fix = atan(-U1_mean_fix / V1);
% 
% %% 2. LOOP DI SENSIBILITÀ (Analysis Mode)
% for i = 1:length(D2_vec)
% 
%     % Geometria Corrente
%     D2_curr = D2_vec(i);
%     R2_curr = D2_curr / 2;
% 
%     % Inizializzazione variabili per il ciclo di convergenza locale
%     eta_loop = 0.85; 
%     rho2_loop = Pt1 * beta_tt / (R * Tt1 * 1.3); % Guess iniziale
% 
%     % Ciclo di convergenza (per trovare b2 e perdite coerenti)
%     % Inizializzazione variabili convergenza locale
%     err_loc = 1;
%     iter_loc = 0;
%     max_iter_loc = 100; % Safety break per non bloccarsi su geometrie impossibili
% 
%     % Ciclo di convergenza (While invece di For)
%     % Inizializzazione variabili convergenza locale
%     err_loc = 1;
%     iter_loc = 0;
%     max_iter_loc = 100; % Safety break per non bloccarsi su geometrie impossibili
% 
%     % Ciclo di convergenza (While invece di For)
%     while abs(err_loc) > 1e-5 && iter_loc < max_iter_loc
%         iter_loc = iter_loc + 1;
% 
% 
%         U2_curr = omega * R2_curr;
% 
% 
%         L_req = cp * Tt1/eta_loop * (beta_tt^((gamma-1)/gamma) - 1);
% 
%         V2_tg_curr = L_req / U2_curr;
%         V2_m_curr = V2_tg_curr / tan(alpha_2);
% 
%         V2_curr = sqrt(V2_m_curr^2 + V2_tg_curr^2);
%         W2_tg_curr = V2_tg_curr - U2_curr;
%         W2_curr = sqrt(V2_m_curr^2 + W2_tg_curr^2);
% 
% 
%         Tt2_curr = Tt1 + L_req/cp;
%         T2_curr = Tt2_curr - V2_curr^2/(2*cp);
%         M2_curr = V2_curr / sqrt(gamma * R * T2_curr);
% 
%         Pt2_curr = Pt1 * (1 + eta_loop * L_req/(cp*Tt1))^(gamma/(gamma-1));
%         P2_curr = Pt2_curr / (1 + (gamma-1)/2 * M2_curr^2)^(gamma/(gamma-1));
% 
%         rho2_new = P2_curr / (R * T2_curr);
%         rho2_loop = 0.5 * rho2_loop + 0.5 * rho2_new; % Rilassamento densità
% 
%         b2_curr = m_dot / (rho2_loop * V2_m_curr * pi * D2_curr);
% 
% 
%         beta2_curr = atan(W2_tg_curr / V2_m_curr);
% 
% 
%         beta_av_curr = 0.5 * (beta1_mean_fix + beta2_curr);
%         N_bl_curr = ceil(2*pi*cos(beta_av_curr)/(0.4*log(D2_curr/(2*R1_m))));
%         if mod(N_bl_curr,2)==1
%             N_bl_curr = N_bl_curr + 1; 
%         end
% 
%         % Spessore pala
%         t_calc = 0.003 * D2_curr;
% 
% 
%         mu_slip = 1 - 0.63*pi/N_bl_curr;
%         V2_tg_inf = (1 - mu_slip)*U2_curr + V2_tg_curr;
%         W2_tg_inf = V2_tg_inf - U2_curr;
%         beta2_geom_curr = atan(W2_tg_inf/V2_m_curr);
% 
%         beta1_geom_mean_curr = atan(1 - t_calc*N_bl_curr/(2*pi*R1_m)*tan(beta1_mean_fix));
% 
%         % losses
%         T_avg_loc = (Tt1 + T2_curr)/2;
%         visc_loc = mu_NH3(T_avg_loc);
% 
%         % 1. Incidence
%         dH_inc = 0.4*(W1_mean_fix - V1/cos(beta1_mean_fix))^2;
% 
%         % 2. Friction
%         D_h_loc = D2_curr*cos(beta2_curr) / (( N_bl_curr/pi + D2_curr*cos(beta2_curr)/b2_curr ) + 0.5 );
%         Re_fr = U2_curr * D_h_loc * rho_t1 / visc_loc;
%         Cf_loc = 0.0412 * Re_fr^(-0.1925);
%         L_hyd = (R2_curr - R1_m)/cos(beta_av_curr); 
%         W_avg_sq = (W1_mean_fix^2 + W2_curr^2)/2;
%         dH_fr = 4 * Cf_loc * W_avg_sq * L_hyd / D_h_loc;
% 
%         % 3. Clearance
%         eps_cl = 5e-4;
%         dH_cl = 0.6 * eps_cl / b2_curr * V2_tg_curr * sqrt(4*pi/(b2_curr*N_bl_curr)*(R1_t_opt^2-R1_h^2)/(R2_curr-R1_t_opt)*V2_tg_curr*V1/(1+rho2_loop/rho1));
% 
%         % 4. Loading & Mixing
%         lb_loc = pi/8 *(D2_curr - R1_m - b2_curr + 0.1)*(2/(cos(beta1_mean_fix)+cos(beta2_curr))); 
%         dW = 2*pi*D2_curr*V2_tg_curr/(N_bl_curr*lb_loc);
%         dH_bl = (dW^2)/48;
% 
%         W_max = 0.5*(W1_mean_fix + W2_curr + dW);
%         D_eq = W_max/W2_curr;
%         if D_eq <= 2
%             W_sep = W2_curr; 
%         else
%             W_sep = W2_curr*D_eq*0.5; 
%         end
%         blockage = 1 - (N_bl_curr*0.002)/(pi*D2_curr);
%         W_out = sqrt((V2_m_curr*blockage)^2 + W2_tg_curr^2);
%         dH_mix = 0.5 * (W_sep - W_out)^2;
% 
%         % 5. Disk Friction
%         visc_disk = mu_NH3(T2_curr);
%         Re_df = U2_curr * R2_curr * rho2_loop / visc_disk;
%         if Re_df < 3e5
%             f_df = 2.67/sqrt(Re_df); 
%         else 
%             f_df = 0.0622/Re_df^0.2; 
%         end
%         dH_disk = f_df * rho2_loop * R2_curr^2 * U2_curr^3 / (8*m_dot);
% 
%         % 6. Recirculation & Leakage
%         D_fac = 1 - W2_curr/W1_mean_fix + 0.75*L_req*W2_curr/((N_bl_curr/pi*(1-2*R1_t_opt/D2_curr)+2*2*R1_t_opt/D2_curr)*W1_mean_fix*U2_curr^2);
%         dH_rec = 0.02 * sqrt(tan(alpha_2)) * D_fac^2 * U2_curr^2;
%         dH_leak = 0.6 * eps_cl/b2_curr * V2_curr * sqrt(4*pi/(b2_curr*N_bl_curr)*((2*R1_t_opt-2*R1_h)/(2*R2_curr-2*R1_t_opt))/(1+rho2_loop/rho1)*V2_tg_curr*V1);
% 
%         % 7. Diffusion & Choke
%         dH_diff = 0; dH_choke = 0;
%         if (W1_tip/W_th > 1.75)
%             dH_diff = 0.5*(W1_tip - 1.75*W_th)^2;
%         end 
% 
% 
%         dH_tot_sens = dH_inc + dH_fr + dH_bl + dH_mix + dH_cl + dH_diff + dH_choke + dH_disk + dH_rec + dH_leak;
% 
%         eta_new = (L_req - dH_tot_sens) / L_req;
% 
% 
%         err_loc = (eta_new - eta_loop)/eta_loop;
% 
%         % Rilassamento
%         eta_loop = eta_loop + 0.5 * (eta_new - eta_loop);
%     end
% 
%     eta_sens(i) = eta_loop;
%     b2_sens(i)  = b2_curr;
%     beta2_sens(i) = beta2_geom_curr * 180/pi;
% end
% 
% %% 3. PLOT RISULTATI
% figure
% plot(D2_vec*1000, eta_sens, 'b-', 'LineWidth', 2)
% hold on
% grid on
% xline(D2*1000, '--k', 'Design Point')
% [max_eta, idx_best] = max(eta_sens)
% plot(D2_vec(idx_best)*1000, max_eta, 'rp', 'MarkerSize', 10, 'MarkerFaceColor','r')
% xlabel('D_2 [mm]')
% ylabel('\eta_rot')
% title('Optimization of Rotor Diameter')
% legend('Sensitivity Curve','Current Design','Optimum Point')
% 
% 
% figure
% plot(D2_vec*1000, b2_sens*1000, 'r-', 'LineWidth', 2)
% hold on
% grid on
% plot(D2_vec(idx_best)*1000, b2_sens(idx_best)*1000, 'ko', 'MarkerFaceColor','k')
% xlabel('D_2')
% ylabel('b_2')
% title('Blade Height Constraint (Mass Conservation)')