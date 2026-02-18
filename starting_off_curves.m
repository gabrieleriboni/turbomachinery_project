close all
clc

% --- 1. Girante (Impeller) ---
% Raggi e Diametri
geom.R1_h = R1_h;           % Raggio Hub Inlet
geom.R1_t = R1_t_opt;       % Raggio Tip Inlet
geom.R1_m = R1_m;           % Raggio Medio Inlet
geom.D2   = D2;             % Diametro Outlet
geom.R2   = D2/2;           % Raggio Outlet

% Larghezze canale (Blade Height)
geom.b1 = b1;               % Altezza pala inlet
geom.b2 = b2;               % Altezza pala outlet

% Caratteristiche Pala
geom.N_bl = N_bl;           % Numero di pale
geom.t    = t;              % Spessore pala
geom.Lz   = Lz;             % Lunghezza assiale (per stima attrito)

geom.beta1_metal_tip  = beta1_geom_tip;   % Angolo metallico tip
geom.beta1_metal_mean = beta1_geom_mean;  % Angolo metallico mean (usato per l'incidenza 1D)
geom.beta1_metal_hub  = beta1_geom_hub;   % Angolo metallico hub
geom.beta2_metal      = beta2_geom;

% --- 2. Diffusore Vaneless ---
geom.D3 = D3;               % Diametro ingresso vaned (fine vaneless)
geom.b3 = b3;               % Larghezza (spesso uguale a b2 o pinchata)

% --- 3. Diffusore Palettato (Vaned Diffuser) ---
geom.N_bl_VD = N_bl_VD;     % Numero pale statore
geom.R3      = R3;          % Raggio ingresso
geom.R4      = R4;          % Raggio uscita (o throat region)
geom.R5      = R5;          % Raggio scarico diffusore
geom.b4      = b4;          % Larghezza

% Angoli metallici statore (Costruzione Wedge/Camber)
geom.alpha3_metal = alpha3_geom; % Angolo ingresso pala statore
geom.alpha4_metal = alpha4_geom; % Angolo uscita pala statore

% Gole e Aree (Critico per il chocking)
geom.W3_throat = W3;             % Larghezza gola ingresso
geom.W4_throat = W4;             % Larghezza gola uscita
geom.Lb_VD     = Lb;             % Lunghezza corda/pala statore
geom.A_th_VD   = N_bl_VD * W3 * b3; % Area di gola fisica totale

% --- 4. Voluta e Cono ---
geom.r_vol_inlet = R5;           % Raggio inizio voluta
geom.D_throat_vol = D_throat;    % Diametro gola voluta (Area A6)
geom.A7_exit      = A7;          % Area uscita cono

% --- 5. Dati di Riferimento ---
geom.omega_design = omega;


%% off design
m_dot_design = 3.0; 
m_dot_vec = linspace(1.5, 4.0, 50); 

% Gas Properties Function (Re-defined for local scope)
mu_NH3 = @(T) mu_0_NH3 * (T./T_0_NH3).^(1.5) .* (T_0_NH3 + S_NH3) ./ (T + S_NH3);

% Pre-allocation
beta_tt_vec = zeros(size(m_dot_vec));
eta_tt_vec  = zeros(size(m_dot_vec));
power_vec   = zeros(size(m_dot_vec));
choke_flag  = zeros(size(m_dot_vec)); % 1 if choked

fprintf('Starting Off-Design Calculation...\n');

for i = 1:length(m_dot_vec)
    m_dot = m_dot_vec(i);

    % 1. IMPELLER INLET (Station 1)
    rho1 = Pt1 / (R * Tt1); % First guess
    err_rho1 = 1;
    iter_in = 0;
    
    while err_rho1 > 1e-5 && iter_in < 50
        % Velocità meridiana (Area geometrica fissa)
        A1_geom = pi * (geom.R1_t^2 - geom.R1_h^2);
        Cm1 = m_dot / (rho1 * A1_geom);
        
        % Temperatura e Pressione Statiche
        T1 = Tt1 - Cm1^2 / (2*cp);
        if T1 < 0, fail_flag=true; break; end % Check fisico
        P1 = Pt1 * (T1/Tt1)^(gamma/(gamma-1));
        
        rho1_new = P1 / (R*T1);
        err_rho1 = abs(rho1_new - rho1)/rho1;
        rho1 = rho1_new;
        iter_in = iter_in + 1;
    end

    % Triangoli di velocità INLET
    U1_tip  = geom.omega_design * geom.R1_t;
    U1_mean = geom.omega_design * geom.R1_m;
    U1_hub  = geom.omega_design * geom.R1_h;
    V1 = Cm1;
    W1_meridional = Cm1;

    W1_tip_t  = -U1_tip;
    W1_mean_t = -U1_mean;
    W1_hub_t  = -U1_hub;

    W1_tip  = sqrt(W1_meridional^2 + W1_tip_t^2);
    W1_mean = sqrt(W1_meridional^2 + W1_mean_t^2);
    W1_hub  = sqrt(W1_meridional^2 + W1_hub_t^2);

    % Angoli di flusso (Flow Angles)
    beta1_flow_tip  = atan(W1_tip_t  / W1_meridional);
    beta1_flow_mean = atan(W1_mean_t / W1_meridional);
    beta1_flow_hub  = atan(W1_hub_t  / W1_meridional);


    % 2. IMPELLER OUTLET (Station 2)