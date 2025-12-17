%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   MAIN CODE, TURBOMACHINERY PROJECT                     %
%                    Ascari Jacopo, Riboni Gabriele                       %
%                        CENTRIFUGAL COMPRESSOR                           %
%                          1: Impeller Inlet                              %
%                          2: Impeller Outlet                             %
%                          3: Vaneless Diffuser                           %
%                          4: Vaned Diffuer                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%% input

beta_tt = 2;
m_dot = 3; %kg/s
Pt1 = 1; %bar
Pt1 = Pt1 * 1e5; % Pa
Tt1 = 300; %K

%%% properties gas %%%
cp_exp = 0.037; % kJ/(mol.K)
cv_exp = 0.028; % kJ/(mol.K)
Mm = 17.03; %g/mol
cp = cp_exp/Mm * 1e6;
cv = cv_exp/Mm * 1e6;
gamma = cp_exp/cv_exp;
R = cp-cv; % kJ/(mol.K)

%%% assumptions %%%

omega = 1600; %rpm with gearbox
% assumiamo inlet con velocità bassa poi iter
% per ora veramente gas perfetto
% v1t =0 axial inlet
alpha_2 = 70;  % piccolo per evitare rircolo come da lezione
alpha_2 = alpha_2 * pi/180;

%% balje

%%% input (omega_s)

dht_is = gamma * R/(gamma -1) * Tt1 * (beta_tt ^((gamma - 1)/gamma) -1);
rho_t1 = Pt1/R/Tt1;
Q_in = m_dot/rho_t1;
omega_s = omega * sqrt(Q_in)/dht_is^(3/4);

%%% output (Ds, psi)

Ds = 4.5;  % first guess  % check su psi per shrouded/un
eta_c = 0.82; %first guess

%% iterative process

%%% input

D2 = Ds * sqrt(Q_in)/(dht_is^(1/4));
L = gamma * R/(gamma -1) * Tt1/eta_c * (beta_tt ^((gamma - 1)/gamma) -1);
psi = L/omega^2/D2^2;  %unshrouded
U2 = omega * D2/2;
V2_tg = L/U2;
tau = V2_tg/U2; %backwards
M2_u = U2/sqrt(gamma * R * Tt1);  % molto vicino al sonic, stare attenti (come es)

rho_vect = [rho_t1];
err_new = 1;
err_rho = [err_new];
i = 0;
iter = [i];
relaxation = 0.2;
while abs(err_new) > 1e-6

    D1_h = 0.2; % assumption
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
    D1_t_opt = R1_t_opt*2;

    T1 = Tt1 - cm1(R1_t_opt).^2/(2*cp);
    P1 = Pt1 * (T1/Tt1)^(gamma/(gamma-1));

    rho_new = P1/R/T1;
    rho_new = rho + relaxation * (rho_new - rho);
    err_new = abs(rho_new - rho)/rho_t1;
    err_rho = [err_rho; err_new];
    rho_vect = [rho_vect; rho_new];
    i = i+1;
    iter = [iter; i];
end
rho1 = rho_vect(end);
%% plot iterative process rho1

figure
semilogy(iter, err_rho, '.-')
grid on
title('err_{\rho}')
xlabel('iter')
ylabel('log(err)')
figure
plot(iter, rho_vect, '.-')
grid on
title('density')
xlabel('iter')
ylabel('err')

%%

nu = D1_h/D1_t_opt;   % check [0.3-0.7]
mer_grad = D1_t_opt/D2; % check [0.5-0.7]
b1 = R1_t_opt - R1_h;
R1_m = 1/2 * (R1_t_opt + R1_h);
D1_m = R1_m * 2;
A1 = pi*D1_m*b1;

% velocity triangles 1
U1_tip = omega * R1_t_opt;
U1_mean = omega * R1_m;
U1_hub = omega * R1_h;

W1_tip_t = -U1_tip;
W1_mean_t = -U1_mean;
W1_hub_t = -U1_hub;
W1_meridional = m_dot/rho/pi/(R1_t_opt^2-R1_h^2);
W1_tip = sqrt(W1_tip_t^2 + W1_meridional^2);
W1_mean = sqrt(W1_mean_t^2 + W1_meridional^2);
W1_hub = sqrt(W1_hub_t^2 + W1_meridional^2);

beta1_tip = atan(W1_tip_t/W1_meridional);
beta1_mean = atan(W1_mean_t/W1_meridional);
beta1_hub = atan(W1_hub_t/W1_meridional);
M1_tip_rel;

V1 = W1_meridional;

% velocity triangles 2

% U2 noto e costante
% V2_tg nota dal lavoro

V2_meridional = V2_tg/tan(alpha_2);
V2 = V2_tg/sin(alpha_2);

W2_tg = V2_tg - U2;
W2_meridional = V2_meridional;
W2 = sqrt(W2_meridional^2 + W2_tg^2);

L_eul = U2*V2_tg;
phi = 4 * Q_in/(pi*U2*D2^2);
%% 2. Impeller Outlet

eta_tt_guess = 0.90;
eta_tt_vect = [eta_tt_guess];
err_eta_new = 1;
err_eta = [err_eta_new];
k = 0;
iter_eta = [k];
relaxation2 = 0.2;

while abs(err_eta_new) > 1e-6
    eta_tt = eta_tt_vect(end);
    DR = abs(W1_mean_t) / W2;

    beta2 = atan(W2_tg/V2_meridional);
    beta2_deg = beta2*180/pi; %torna unshrouded

    Tt2 = L/cp + Tt1;
    T2 = Tt2 - V2^2/(2*cp);
    M2 = V2/sqrt(gamma * R *T2);

    chi = cp*(T2-T1)/L;



    Pt2 = Pt1 * (1 + (eta_tt*L)/(cp*Tt1))^(gamma/(gamma - 1));
    P2 = Pt2/(1 + (gamma-1)/2 * M2^2)^(gamma/(gamma - 1));
    rho2 = P2/R/T2;

    b2 = m_dot/(rho2*V2_meridional*pi*D2);
    AR = b2/D2; % check [0.03-0.15]

    solidity = 1/0.4; % c/s, da eckert suggest
    dtheta = log(D2/D1_m) * tan(beta2);
    R2 = D2/2;
    c = R2 * dtheta/(sin(beta2));
    s = c/solidity;

    beta_av = 0.5 * (beta1_mean + beta2);

    % Stodola formula
    N_bl = ceil(2*pi*cos(beta_av)/(0.4*log(D2/D1_m)));

    if mod(N_bl,2)==1
        N_bl = N_bl + 1;
    end


    %mu = 1 - sqrt(cos(beta2))/N_bl^0.7;
    mu = 1 - 0.63*pi/N_bl; % stanitz
    V2_tg_inf = (1- mu)*U2 + V2_tg;
    W2_tg_inf = V2_tg_inf - U2;
    beta2_geom = atan(W2_tg_inf/W2_meridional);
    beta2_geom_deg = beta2_geom *180/pi; % va bene per la correlazione

    t = 3e-3;
    % For ammonia applications, avoid going
    % below 2.5 mm – 3.0 mm at the tip.
    % Standard aero-compressors might go down to 0.8 mm,
    % but this is too fragile for industrial ammonia service where
    % impurities or liquid droplets might exist
    beta1_geom_tip = atan(1-t*N_bl/(2*pi*R1_t_opt)*tan(beta1_tip));
    beta1_geom_mean = atan(1-t*N_bl/(2*pi*R1_m)*tan(beta1_mean));
    beta1_geom_hub = atan(1-t*N_bl/(2*pi*R1_h)*tan(beta1_hub));

    %% rotor losses

    % incidence=0 per costruzione ma non off-design
    
    %%% IMPELLER INTERNAL
    % skin friction (Jansen, 1967)
    visc_din1 = 1.76e-05; %sutherland
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
    % Aungier:
    % m_dot_cl = ;
    % D_avg = 0.5*(D2+ D1_m);
    % b_avg = 0.5*(b1+b2);
    % l_meridional = D2/2 - (D_e);
    % dP_cl = m_dot*(D2*V2_tg)/(N_bl*D_avg*b_avg*l_meridional);
    % dH_cl = m_dot_cl*dP_cl/(m_dot*rho1);

    % blade loading (Aungier, 1995) assumption of no inlet swirl (come da
    % china)
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
    pitch = pi*D1_m/N_bl - t;
    A_th_geom = pitch * cos(beta1_geom_mean)*b1;
    A_th = 0.97 * A_th_geom;
    W_th = m_dot/rho1/A_th;
    dH_diff = 0.4*(W1_mean - W_th)^2;
    if (W1_tip/W_th > 1.75) && (dH_diff < 0.5*(W1_tip - 1.75*W_th)^2)
            dH_diff = 0.5*(W1_tip - 1.75*W_th)^2;
    end
   
    % % choke losses
    
    % A_th_star = m_dot/rho1/sqrt(gamma*R*T1);
    % C_r = sqrt(A1*cos(beta1_geom_mean)/A_th);
    % if C_r > 1-(A1*cos(beta1_geom_mean)/A_th -1)^2
    %     C_r = 1-(A1*cos(beta1_geom_mean)/A_th -1)^2;
    % end
    % X = 11 - 10*(C_r*A_th)/A_th_star;
    % if X<=0
    %     dH_choke = 0;
    % else 
    %     dH_choke = 1/2*W1_mean^2*(0.05*X+X^7);
    % end

    

    %%% PARASSITIC
    % disk friction  (Daily and Nece)
    visc_din2 = 2.00e-05; %sutherland
    Re_df = U2*R2/visc_din2*rho2;
    if Re_df < 3e5
        f_df = 2.67/Re_df^0.5;
    else
        f_df = 0.0622/Re_df^0.2;
    end
    dH_disk = f_df * (rho1 + rho2)*R2^2*U2^3/(8*m_dot);
    %dH_disk = f_df * (1+P2/P1) * 0.5 * L_eul * V1/U2 * R2^2/R1_t_opt*(1-(R1_h/R2)^2);

    % recirculation (Coppage)
    D_factor = 1 - W2/W1_tip + 0.75*L_eul*W2/((N_bl/pi*(1-D1_t_opt/D2)+2*D1_t_opt/D2)*W1_tip*U2^2);
    dH_rec = 0.02*sqrt(tan(alpha_2))*D_factor^2*U2^2;

    % leakage (Jansen)
    dH_leak = 0.6 * eps/b2*V2*sqrt(4*pi/(b2*N_bl)*((R1_t_opt-R1_h)/(R2-R1_t_opt))/(1+rho2/rho)*V2_tg*V1);

    

    % AGGIUNGERE dH_diff
    dH_tot_internal =  dH_fr + + dH_bl + dH_mix + dH_cl + dH_diff;
    dH_tot_parassitic = dH_disk + dH_rec + dH_leak;
    dH_tot = dH_tot_internal + dH_tot_parassitic;

    dH_t_new = L - dH_tot;

    eta_tt_new = dH_t_new/L;

    eta_tt_new = eta_tt + relaxation2 * (eta_tt_new - eta_tt);
    eta_tt_vect = [eta_tt_vect; eta_tt_new];

    err_eta_new = abs(eta_tt_new - eta_tt)/eta_tt;
    err_eta = [err_eta; err_eta_new];
 
    k = k+1;
    iter_eta = [iter_eta; k];
end

%% plot iterative process eta_tt

figure
semilogy(iter_eta, err_eta, '.-')
grid on
title('err_{\eta}')
xlabel('iter')
ylabel('log(err)')
figure
plot(iter_eta, eta_tt_vect, '.-')
grid on
title('\eta')
xlabel('iter')
ylabel('err')
% itero per diff di eta poi ricalcolo Pt2_new e P2_new e b2_new

%% 3. Vaneless Diffuser

%b3 = % we introduce a slight pinch
D3 = 1.2*D2;   % NOTARE CHE D5/D2 = 1.5/1.6 MINIMUM
R3 = D3/2;
Tt3 = Tt2;
visc_din3 = visc_din2;


% rho3 = rho2;  %first guess
% V3 = V2;      %first guess
% k = 0.01;
% rho_avg = 0.5*(rho2 + rho3);
% V_avg = 0.5 * (V2 + V3);
% D_avg = 0.5 * (D2+D3);
% Re_avg = rho_avg * V_avg * D_avg/visc_din3;
% Cf_vaned = k*(1.8e5/Re_avg)^0.2;
% V3_tg = V2_tg/(D3/D2+0.5*pi*Cf_vaned*rho2*V2_tg*D3*(D3-D2)/m_dot);
% V3_meridional = (V2_meridional*rho2*pi*D2*b2)/(rho3*pi*D3*b3);
% V3_new = sqrt(V3_meridional^2+V3_tg^2);

% dH_t_VD = Cf_vaned * V1^2 * R2 *(1-(R2/R3)^1.5)/(1.5*b2*cos(alpha_2));
% 
% T3 = Tt3-1/(2*cp)*V3^2;
% M3 = V3/sqrt(gamma*R*T3);
% Tt3_is = Tt3-dH_t_VD/cp;
% Pt3 = Pt2*(Tt3_is/Tt2)^(gamma/(gamma-1));
% P3 = Pt3/(1+(gamma-1)/2*M3^2)^(gamma/(gamma-1));
% rho3_new = P3/(T3*R);

% fare ciclo while iterando su rho, poi fare vaned diffuser

%eta_d = 0.9;  %first guess
% P4 = P3 + Cp*(Pt3-P3);
% T4_is = T3*(P4/P3)^((gamma-1)/gamma);
% T4 = (T4_is-T3)/eta_d_guess+T3;
% Tt4 = Tt2;
% V4 = sqrt(2*cp*(Tt4-T4));
% M4 = V4/sqrt(gamma*R*T4);
% rho4 = P4/(R*T4);

% ee = ;
% AR = ;
% N_bl_d = ceil(360/ee);
% 
% H3 = (2*pi*R3*cos(alpha3))/N_bl_d;
% S3 = H3*b3;
% S4 = AR*S3;

% L3_b = ;
% Ld = b3*L3_b;
% alpha4 = atan(R3*sin(alpha3)/(Ld+R3*cos(alpha3)));
% R4 = sqrt(Ld^2+R3^2+2*Ld*R3*cos(alpha3));
% H4 = (2*pi*R4*cos(alpha4))/N_bl_d;
% b4 = S4/H4;

% Pt4_new = P4*(1+(gamma-1)/2*M4^2)^(gamma/(gamma-1));
% L_is_new = Cp*T_t1*((Pt4_new/Pt1)^((gamma-1)/gamma)-1);
% 
% eta_new = L_is_new/L;

% fare ciclo while su tutto eta


%% Checks

m_dot1 = rho * V1 * pi * (R1_t_opt^2 - R1_h^2); % V1_meridional = V1
m_dot2 = rho2 * V2_meridional * D2 * b2 * pi;
%m_dot3 = rho3 * V3_meridional*D3*b3*pi;
%m_dot4 = rho4*V4*S4*N_bl_d;

check = {'psi', psi * 4,'[0.7-0.8] unshrouded ([0.6-0.7] shrouded)';
    'tau', tau, '<1 backwards';
    'M2_u', M2_u, '<1 subsonic'
    'nu', nu, '[0.3-0.7]';
    'mer_grad', mer_grad, '[0.5-0.7]';
    'M1_tip_rel',M1_tip_rel,'<1.3';
    'DR',DR, '[]';
    'M2',M2, '[<1.3]';
    'chi',chi, '[0.5-0.6] unshrouded ([0.6-0.7] shrouded)';
    'beta2', beta2_geom_deg, '[-20 -30] unshrouded ([-40 -50] shrouded)';
    'AR', AR, '[0.03-0.15]';
    'blade height (nu)', nu, '[0.3-0.7]';
    'alpha2', rad2deg(alpha_2), '[<80 avoiding recirculation]';
    'mer_grad', mer_grad, '[0.5-0.7]';
    'U2', U2, '[for mechanical stresses]';
    'm. flow in',m_dot,'3';
    'm. flow 1',m_dot1,'3';
    'm. flow 2',m_dot2,'3';
    'm. flow 3','m_dot3','3';
    'm. flow 4','m_dot4','3'}

%% plot velocity triangles

figure
subplot(1,2,1)
title('impeller inlet')
xlabel('meridional')
ylabel('tangential')
hold on

% W1 parte da U1
quiver(0,U1_tip,  W1_meridional,W1_tip_t, 0, 'r')
quiver(0,U1_mean, W1_meridional,W1_mean_t, 0, 'r')
quiver(0,U1_hub,  W1_meridional,W1_hub_t, 0, 'r')

% U1
quiver(0,0, 0,U1_tip, 0,  'k')
quiver(0,0, 0,U1_mean, 0, 'k')
quiver(0,0, 0,U1_hub, 0, 'k')

% V1
quiver(0,0, V1,0, 0,'c')
axis equal
grid on


offset = 1.5;
text(V1/2, offset, 'V_1', 'Color','c');
text(0, U1_tip/2 + offset, 'U_{1,tip}', 'Color','k');
text(0, U1_mean/2 + offset, 'U_{1,mean}', 'Color','k');
text(0, U1_hub/2 + offset, 'U_{1,hub}', 'Color','k');
text(W1_meridional/2, U1_tip + W1_tip_t/2 + offset, 'W_{1,tip}', 'Color','r');
text(W1_meridional/2, U1_mean + W1_mean_t/2 + offset, 'W_{1,mean}', 'Color','r');
text(W1_meridional/2, U1_hub + W1_hub_t/2 + offset, 'W_{1,hub}', 'Color','r');

subplot(1,2,2)
title('impeller outlet')
xlabel('meridional')
ylabel('tangential')
hold on

% W2 parte da U2
quiver(0,U2,  W2_meridional,W2_tg, 0, 'r')

% U2
quiver(0,0, 0,U2, 0, 'k')

% V2
quiver(0,0, V2_meridional, V2_tg, 0,'c')
axis equal
grid on

offset2 = 1.5;
text(V2_meridional/2 + offset2, V2_tg/2 + offset2, 'V_2', 'Color','c');
text(0, U2/2 + offset2, 'U_2', 'Color','k');
text(W2_meridional/2, U2 + W2_tg/2 + offset2, 'W_2', 'Color','r');


