%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     COMPRESSOR GEOMETRY GENERATOR                       %
%                 2D Meridional Profile Visualization                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

% 1. Run the main design to get the workspace variables
fprintf('Running main design script to extract geometry...\n');
evalc('main_def'); % Runs main_def without printing its output to console
fprintf('Design complete. Extracting geometric parameters...\n\n');

%% --- Extract Geometric Points ---

% Convert all dimensions to millimeters for better plotting
scale = 1000; 

% 1. Impeller Inlet (Station 1)
r1_hub  = R1_h * scale;
r1_tip  = R1_t_opt * scale;
b1_val  = b1 * scale;

% 2. Impeller Outlet (Station 2)
r2      = R2 * scale;
b2_val  = b2 * scale;

% 3. Vaneless Diffuser Outlet (Station 3)
r3      = R3 * scale;
b3_val  = b3 * scale;

% 4. Vaned Diffuser Outlet (Station 4)
r4      = R4 * scale;
b4_val  = b4 * scale;

% 5. Second Vaneless Diffuser Outlet (Station 5)
r5      = R5 * scale;
b5_val  = b5 * scale;

% 6. Volute Mean Radius (Station 6)
r6_mean = r6 * scale;
r_volute_cross = r_sect * scale; % Radius of the volute circular cross-section

%% --- Estimate Axial Dimensions (Z-axis) ---
% Since this is a 1D code, we use standard empirical proportions to estimate 
% the axial length of the impeller for visualization.
L_ax = 0.3 * (r2 * 2); % Typical axial length ~ 30% of D2

% Z-Coordinates for the Shroud (Outer casing)
z_inlet = 0;
z_imp_out = L_ax;

%% --- Construct the Meridional Curves ---

% --- HUB CURVE ---
% Inlet hub point
hub_z(1) = z_inlet;
hub_r(1) = r1_hub;
% Impeller outlet hub point
hub_z(2) = z_imp_out;
hub_r(2) = r2;
% Diffusers hub line (straight radial)
hub_z(3) = z_imp_out;
hub_r(3) = r5;

% Generate a smooth curve for the Impeller Hub using a spline
% Fix: ensure strictly increasing x-values for the spline control points
z_mid = hub_z(1) + 0.3 * (hub_z(2) - hub_z(1));
z_hub_spline = linspace(hub_z(1), hub_z(2), 50);
r_hub_spline = spline([hub_z(1), z_mid, hub_z(2)], [hub_r(1), hub_r(1), hub_r(2)], z_hub_spline);

% --- SHROUD CURVE ---
% Inlet tip point
shroud_z(1) = z_inlet;
shroud_r(1) = r1_tip;
% Impeller outlet tip point (shifted by b2)
shroud_z(2) = z_imp_out - b2_val;
shroud_r(2) = r2;
% Vaneless Diffuser tip (pinch applied)
shroud_z(3) = z_imp_out - b3_val;
shroud_r(3) = r3;
% Vaned Diffuser tip (constant width)
shroud_z(4) = z_imp_out - b4_val;
shroud_r(4) = r4;
% Second Vaneless tip
shroud_z(5) = z_imp_out - b5_val;
shroud_r(5) = r5;

% Generate a smooth curve for the Impeller Shroud using a spline
z_mid_s = shroud_z(1) + 0.6 * (shroud_z(2) - shroud_z(1));
z_shroud_spline = linspace(shroud_z(1), shroud_z(2), 50);
r_shroud_spline = spline([shroud_z(1), z_mid_s, shroud_z(2)], [shroud_r(1), shroud_r(1)*1.2, shroud_r(2)], z_shroud_spline);

%% --- Plotting the Meridional Profile ---

figure('Name', 'Compressor Meridional Geometry', 'Position', [100, 100, 900, 600])
hold on; grid on; axis equal;
title('\textbf{Centrifugal Compressor - Meridional Profile}', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Axial Direction $Z$ [mm]', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('Radial Direction $R$ [mm]', 'Interpreter', 'latex', 'FontSize', 12)

% 1. Plot the Impeller
plot(z_hub_spline, r_hub_spline, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Hub Curve');
plot(z_shroud_spline, r_shroud_spline, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Shroud Curve');

% Draw Inlet Line
plot([hub_z(1), shroud_z(1)], [hub_r(1), shroud_r(1)], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Inlet (St. 1)');

% Draw Impeller Outlet Line
plot([hub_z(2), shroud_z(2)], [hub_r(2), shroud_r(2)], 'g--', 'LineWidth', 1.5, 'DisplayName', 'Impeller Exit (St. 2)');

% 2. Plot the Diffuser Sections (Hub side is flat, Shroud side steps/tapers)
% Vaneless Diffuser 1
plot([hub_z(2), hub_z(3)], [r2, r3], 'k-', 'LineWidth', 2.5, 'HandleVisibility','off');
plot([shroud_z(2), shroud_z(3)], [r2, r3], 'b-', 'LineWidth', 2.5, 'HandleVisibility','off');
plot([hub_z(3), shroud_z(3)], [r3, r3], 'm:', 'LineWidth', 1.5, 'DisplayName', 'Vaneless Exit (St. 3)');

% Vaned Diffuser
plot([hub_z(3), hub_z(3)], [r3, r4], 'k-', 'LineWidth', 2.5, 'HandleVisibility','off');
plot([shroud_z(3), shroud_z(4)], [r3, r4], 'b-', 'LineWidth', 2.5, 'HandleVisibility','off');
% Draw Vaned Diffuser Domain (Shaded region to represent blades)
patch([shroud_z(3), hub_z(3), hub_z(3), shroud_z(4)], [r3, r3, r4, r4], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Vaned Diffuser Domain');
plot([hub_z(3), shroud_z(4)], [r4, r4], 'c:', 'LineWidth', 1.5, 'DisplayName', 'Vaned Diff. Exit (St. 4)');

% Second Vaneless Diffuser
plot([hub_z(3), hub_z(3)], [r4, r5], 'k-', 'LineWidth', 2.5, 'HandleVisibility','off');
plot([shroud_z(4), shroud_z(5)], [r4, r5], 'b-', 'LineWidth', 2.5, 'HandleVisibility','off');
plot([hub_z(3), shroud_z(5)], [r5, r5], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Diffuser Exit (St. 5)');

% 3. Plot the Volute Cross Section
% Assume the volute is a circle centered at r6_mean
theta_circle = linspace(0, 2*pi, 100);
volute_z = z_imp_out - b5_val/2 + r_volute_cross * cos(theta_circle); % Center axially roughly with the diffuser exit
volute_r = r6_mean + r_volute_cross * sin(theta_circle);
plot(volute_z, volute_r, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, 'DisplayName', 'Volute Cross-Section (St. 6)');

% Centerline (Axis of Rotation)
plot([-50, max(volute_z)+50], [0, 0], 'k-.', 'LineWidth', 1, 'DisplayName', 'Axis of Rotation');

% Formatting
legend('Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 10);
ylim([0, max(volute_r) * 1.1]);
xlim([-50, max(volute_z) * 1.5]);

% --- Print Geometric Summary ---
fprintf('==========================================\n');
fprintf('        GEOMETRIC SUMMARY (mm)            \n');
fprintf('==========================================\n');
fprintf('--- Impeller ---\n');
fprintf('Inlet Hub Dia (D1_h):  %6.2f\n', r1_hub*2);
fprintf('Inlet Tip Dia (D1_t):  %6.2f\n', r1_tip*2);
fprintf('Outlet Dia    (D2):    %6.2f\n', r2*2);
fprintf('Outlet Width  (b2):    %6.2f\n', b2_val);
fprintf('Blade Count   (N_bl):  %d\n', N_bl);
fprintf('--- Diffusers ---\n');
fprintf('Vaneless Out Dia (D3): %6.2f\n', r3*2);
fprintf('Vaned Diff Out Dia(D4):%6.2f\n', r4*2);
fprintf('Vaned Diff Blades:     %d\n', N_bl_VD);
fprintf('Diffuser Exit Dia(D5): %6.2f\n', r5*2);
fprintf('Exit Width    (b5):    %6.2f\n', b5_val);
fprintf('--- Volute ---\n');
fprintf('Mean Radius   (r6):    %6.2f\n', r6_mean);
fprintf('Section Dia   (D_th):  %6.2f\n', r_volute_cross*2);
fprintf('==========================================\n');
