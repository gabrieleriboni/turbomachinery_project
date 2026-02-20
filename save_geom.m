%% save_geometry_from_design.m
% Run this AFTER your MAIN design script has converged (variables in workspace)

% --- Pack geometry (fixed) + design-point data into structs ---
geom = struct();

% --- Inlet / impeller inlet geometry (fixed) ---
geom.D1_h      = D1_h;
geom.R1_h      = R1_h;
geom.R1_t      = R1_t_opt;
geom.D1_t      = D1_t_opt;
geom.R1_m      = R1_m;
geom.D1_m      = D1_m;
geom.b1        = b1;
geom.A1        = A1;

% --- Impeller outlet geometry (fixed) ---
geom.D2        = D2;
geom.R2        = R2;
geom.b2        = b2;
geom.U2        = U2;
geom.N_bl      = N_bl;
geom.t         = t;
geom.t_te      = t_te;
geom.eps       = eps;

% NOTE: beta2_geom is computed in your loop and is the one you want to freeze
geom.beta2_geom = beta2_geom;

% --- Impeller inlet geometric angles (fixed) ---
geom.beta1_geom_tip  = beta1_geom_tip;
geom.beta1_geom_mean = beta1_geom_mean;
geom.beta1_geom_hub  = beta1_geom_hub;
geom.A_th_imp_geom   = A_th_geom;

% --- Vaneless diffuser geometry (fixed) ---
geom.R3        = R3;
geom.D3        = D3;
geom.b3        = b3;

% --- Vaned diffuser geometry (fixed) ---
geom.N_bl_VD     = N_bl_VD;
geom.R4          = R4;
geom.D4          = D4;
geom.b4          = b4;
geom.alpha3_geom = alpha3_geom;
geom.alpha4_geom = alpha4_geom;
geom.W3          = W3;
geom.W4          = W4;
geom.Lb          = Lb;
geom.theta_c     = theta_c;
geom.camber      = camber;
geom.sigma       = sigma;
geom.theta       = theta;

% --- Second vaneless diffuser geometry (fixed) ---
geom.R5        = R5;
geom.D5        = D5;
geom.b5        = b5;

% --- Volute geometry (we freeze design throat area + mean radius) ---
geom.SP         = SP;
geom.K_friction = K_friction;
geom.r6         = r6;         % mean radius (fixed)
geom.A6         = A6;         % throat area (fixed)
geom.D_throat   = D_throat;   % derived at design (fixed)

% --- Exit cone geometry (fixed) ---
geom.AR_cone = AR_cone;

% --- Operating condition & gas props (store so off-design script is self-consistent) ---
geom.m_dot_design  = m_dot;
geom.beta_tt_design= beta_tt;
geom.Pt1           = Pt1;
geom.Tt1           = Tt1;

geom.cp_exp = cp_exp;
geom.cv_exp = cv_exp;
geom.Mm     = Mm;
geom.cp     = cp;
geom.cv     = cv;
geom.gamma  = gamma;
geom.Rgas   = R;

geom.rpm   = rpm;
geom.omega = omega;

% --- Design-point outputs (for plot anchoring) ---
design = struct();
design.err_final     = err_tot_vect(end);
design.iter_final    = iter_tot(end);
design.eta_c         = eta_c_vect(end);
design.eta_stage_tt  = eta_stage_tt;
design.Pt7           = Pt7;
design.beta_tt       = Pt7/Pt1;  % should match beta_tt target (≈2)

save('geom_design.mat','geom','design');

fprintf('\nSaved geom_design.mat\n');
fprintf('Design point: m_dot=%.6f kg/s | beta_tt=%.6f | eta_c=%.6f | eta_stage=%.6f\n', ...
    geom.m_dot_design, design.beta_tt, design.eta_c, design.eta_stage_tt);
