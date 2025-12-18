%% EGB322 Practical Calculations
close all; clear; clc;

% --- To Do ---
% 1. Update states summary to include temp.

% --- Configuration and Data Loading ---
% Read csv file
engine_data_table = readtable('lab_table_filled.csv');
% To initialise things
num_rows = size(engine_data_table, 1);

% Load data from the Excel file, only getting first 2 columns
data_1400RPM_diesel = readtable('Diesel_PV_Data_1400RPM.xlsx', 'Range', 'A:B');
data_2000RPM_diesel = readtable('Diesel_PV_Data_2000RPM.xlsx', 'Range', 'A:B');
data_1400RPM_petrol = readtable('Petrol_PV_Data_1400RPM.xlsx', 'Range', 'A:B');

% Extract data
V_diesel_1400 = data_1400RPM_diesel.Volume;
P_diesel_1400 = data_1400RPM_diesel.Pressure;
V_diesel_2000 = data_2000RPM_diesel.Volume;
P_diesel_2000 = data_2000RPM_diesel.Pressure;
V_petrol_1400 = data_1400RPM_petrol.Volume;
P_petrol_1400 = data_1400RPM_petrol.Pressure;

% Rename engine_data_table columns for easy dot-notation access
engine_data_table.Properties.VariableNames{'NominalSpeed_RPM_'}  = 'NominalSpeed';
engine_data_table.Properties.VariableNames{'EngineLoad___'}      = 'EngineLoad';
engine_data_table.Properties.VariableNames{'MeasuredSpeed_RPM_'} = 'MeasuredSpeed';
engine_data_table.Properties.VariableNames{'Throttle___'}        = 'Throttle';
engine_data_table.Properties.VariableNames{'Torque_Nm_'}         = 'Torque';
engine_data_table.Properties.VariableNames{'Power_kW_'}          = 'Power';
engine_data_table.Properties.VariableNames{'ElevationDiff_mm_'}  = 'ElevationDiff';
engine_data_table.Properties.VariableNames{'FuelRate_LPM_'}      = 'FuelRate';
engine_data_table.Properties.VariableNames{'CoolRate_LPM_'}      = 'CoolRate';
engine_data_table.Properties.VariableNames{'WaterIn__C_'}        = 'WaterIn';
engine_data_table.Properties.VariableNames{'WaterOut__C_'}       = 'WaterOut';
engine_data_table.Properties.VariableNames{'CalIn__C_'}          = 'CalIn';
engine_data_table.Properties.VariableNames{'CalOut__C_'}         = 'CalOut';
engine_data_table.Properties.VariableNames{'ExhIn__C_'}          = 'ExhIn';
engine_data_table.Properties.VariableNames{'ExhOut__C_'}         = 'ExhOut';

% Display the tables
disp(' ')
disp('Raw Data Table')
disp(engine_data_table)
disp(' ')

% --- Figure Flags ---
% Set flags to draw and save figures
draw_engine_performance_plots = true;
draw_heat_balance_plots = true;
draw_PV_diagrams_plots = true;
draw_diesel_petrol_plot = true;
save_plots = false;

% --- Engine Specifications and Parameters ---
% Engine specifications
bore_diameter = 0.084; % 84 mm [m]
bore_length = 0.1;     % 100 mm [m]
number_cylinders = 4;  % Number of engine cylinders
n_r_strokes = 2; % Number of crankshaft cycles per power cycle (4-stroke)
% Volumetric displacement of engine [m^3]
V_d = (pi * bore_diameter^2 / 4) * bore_length * number_cylinders;

% Fuel and Fluid parameters
density_fuel = 860;  % kg m^-3 (Diesel)
Q_HV_diesel = 44000; % kJ kg^-1 (Heating Value)
density_air = 1.184; % kg m^-3 (Density at 25 C)
c_p_water = 4.2;     % kJ kg^-1 K^-1

% Air/Gas specific parameters for PV analysis
c_p_air = 1.005;     % kJ kg^-1 K^-1, at 300 K
k_air = 1.4;         % Isentropic index for air
T1 = 24.7 + 273.15;  % K, Assumed Temperature at State 1

%% Engine Performace and Heat Balance
% Specific heat for exhaust gas assumed same as air at average exhasut
% tempreture
c_p_exhaust_vector = zeros(num_rows, 1);
c_p_exhaust_vector(2) = 1.022; % Row 2 specific
c_p_exhaust_vector(6) = 1.028; % Row 6 specific

% Initialize storage for results
results_engine_performance = zeros(num_rows, 14);
% Heat Balance is only calculated for rows 2 and 6
heat_balance_rows = [2, 6]; % For indexing
results_heat_balance = [];

% Loop through rows in imported data to do calculations
for i = 1:num_rows
    row = engine_data_table(i, :); 
    c_p_exhaust = c_p_exhaust_vector(i);
    
    [performance_row, heat_balance_row] = ...
        calculate_data_row(row, V_d, density_fuel, density_air, ...
        Q_HV_diesel, c_p_water, c_p_exhaust, n_r_strokes);
    
    results_engine_performance(i, :) = performance_row;
    
    if ismember(i, heat_balance_rows) % Only do on set rows
        results_heat_balance = [results_heat_balance; heat_balance_row];
    end
end

% --- Display Table of Results ---
% Engine Performance Table
column_names_engine_performance = {
    'Nominal_RPM', 'Engine_Load_pct', 'Measured_RPM', 'Power_kW', ...
    'Torque_Nm', 'BMEP_MPa', 'SFC_g_kW_hr', 'Air_Fuel_Ratio', ...
    'Vol_Efficiency_pct', 'Thermal_Efficiency_pct', 'Vol_Flow_Fuel_m3_s', ...
    'Mass_Flow_Fuel_kg_s', 'Vol_Flow_Air_m3_s', 'Mass_Flow_Air_kg_s'};
engine_performance_table = array2table(results_engine_performance, ...
    'VariableNames', column_names_engine_performance);
disp(' ');
disp('Engine Performance Summary Table');
disp(engine_performance_table);
disp(' ');

% Heat Balance Table
column_names_heat_balance = {
    'Nominal_RPM', 'Engine_Load_pct', 'Q_fuel_kW', 'Q_cooling_kW', 'Q_exhaust_kW', ...
    'Brake_Power_kW', 'Q_other_kW', 'Mass_Flow_Exhaust_kg_s', 'Mass_Flow_Water_kg_s'};
heat_balance_table = array2table(results_heat_balance, ...
    'VariableNames', column_names_heat_balance);
disp(' ');
disp('Heat Balance Summary Table');
disp(heat_balance_table);
disp(' ');

% --- Plotting ---
if draw_engine_performance_plots
    plot_engine_performance(engine_performance_table, 'figures/engine_performance_1', save_plots);
end
if draw_heat_balance_plots
    plot_heat_balance_pie_chart(1, heat_balance_table, 'figures/energy_balance_1', save_plots); % Row 2 data
    plot_heat_balance_pie_chart(2, heat_balance_table, 'figures/energy_balance_2', save_plots); % Row 6 data
end

%% PV Analysis
% To get the polytropic coefficients, visually fit the curves on the
% log-log plots using the fitted_points vectors below. 
%
% LHS Points: Pick anywhere on the linear section for both, this
% isn't critical.
%
% RHS points: Have both with the same volume. Pick the largest volume that
% is on both linear sections, as V1_comp will be used for V_max = V1 = V4.

% --- 1400 RPM ---
row_index_1400 = 2;
fitted_points_1400_diesel = [
    111, 4.7e-4, ...  % P1_comp, V1_comp, Bottom Right
    2678, 4e-5, ...   % P2_comp, V2_comp, Bottom Left
    2524, 1e-4, ...   % P1_exp,  V1_exp,  Top Left
    386, 4.7e-4       % P2_exp,  V2_exp,  Top Right
];

[PV_analysis_1400RPM, states_1400RPM] = process_pv_data(V_diesel_1400, P_diesel_1400, fitted_points_1400_diesel, ...
    engine_performance_table, row_index_1400, n_r_strokes, number_cylinders, V_d, ...
    Q_HV_diesel, c_p_air, k_air, T1, draw_PV_diagrams_plots, save_plots);

% --- 2000 RPM ---
row_index_2000 = 6;
fitted_points_2000_diesel = [
    145.6, 5.e-4, ... % P1_comp, V1_comp, Bottom Right
    2209, 7e-5, ...   % P2_comp, V2_comp, Bottom Left
    3101, 1e-4, ...   % P1_exp,  V1_exp,  Top Left
    403, 5e-4         % P2_exp,  V2_exp,  Top Right
];

[PV_analysis_2000RPM, states_2000RPM] = process_pv_data(V_diesel_2000, P_diesel_2000, fitted_points_2000_diesel, ...
    engine_performance_table, row_index_2000, n_r_strokes, number_cylinders, V_d, ...
    Q_HV_diesel, c_p_air, k_air, T1, draw_PV_diagrams_plots, save_plots);

% Combine the results for both RPMs into a single table
combined_PV_analysis = [PV_analysis_1400RPM; PV_analysis_2000RPM];
column_names_PV = {
    'Nominal_RPM', 'Measured_RPM', 'Indicated_Work_kJ', 'Indicated_Power_kW', ...
    'Mechanical_Efficiency_pct', 'IMEP_MPa', 'Heat_Input_Fitted_pct', ...
    'Cut_off_Ratio', 'Energy_In_kJ', 'Energy_Fuel_kJ', 'Mass_Fuel_kg', 'Mass_Air_kg'
    };
PV_analysis_table_combined = array2table(combined_PV_analysis, 'VariableNames', column_names_PV);

disp(' ');
disp('PV Analysis Summary Table');
disp(PV_analysis_table_combined);
disp(' ');

% PV States Table
column_names_states = {'Nominal_RPM', 'n_comp', 'C_comp', 'n_exp', 'C_exp', 'P1_kPa', 'V1_m3', 'T1_k', 'P2_kPa', 'V2_m3', 'T2_k', 'P3_kPa', 'V3_m3', 'T3_k', 'P4_kPa', 'V4_m3', 'T4_k'};
states_combined = [states_1400RPM; states_2000RPM];
PV_states_table = array2table(states_combined, 'VariableNames', column_names_states);

disp(' ');
disp('PV States Summary Table');
disp(PV_states_table);
disp(' ');


%% Diesel and Petrol PV Plots at 1400 RPM
% Uses fitted_points_1400_diesel from above for the diesel section.

% Set this flag to display the log-log plot for fitting the curve to get
% the polytropic coefficients for the petrol cycle.
draw_petrol_log_log_plot = false;

fitted_points_1400_petrol = [
    116, 4.28e-4, ...   % P1_comp, V1_comp, Bottom Right
    973, 6.73e-5, ...   % P2_comp, V2_comp, Bottom Left
    1438, 1.057e-4, ... % P1_exp,  V1_exp,  Top Left
    222, 4.3e-4         % P2_exp,  V2_exp,  Top Right
];

[n_comp_petrol, C_comp_petrol, n_exp_petrol, C_exp_petrol, states_petrol] = polytropic_coefficient_fitting(P_petrol_1400, V_petrol_1400, fitted_points_1400_petrol, 1400, T1, k_air, draw_petrol_log_log_plot, 'petrol_1400RPM', save_plots);
%plot_fitted_PV_curve(P_petrol_1400, V_petrol_1400, n_comp_petrol, C_comp_petrol, n_exp_petrol, C_exp_petrol, states_petrol, 1400, 'petrol_1400RPM', save_plots);

[n_comp_diesel, C_comp_diesel, n_exp_diesel, C_exp_diesel, states_diesel] = polytropic_coefficient_fitting(P_diesel_1400, V_diesel_1400, fitted_points_1400_diesel, 1400, T1, k_air, false, 'diesel_1400RPM', save_plots);
%plot_fitted_PV_curve(P_diesel_1400, V_diesel_1400, n_comp_diesel, C_comp_diesel, n_exp_diesel, C_exp_diesel, states_diesel, 1400, 'diesel_1400RPM', save_plots);

if draw_diesel_petrol_plot
    plot_fitted_PV_curve_dual_engine(P_diesel_1400, V_diesel_1400, states_diesel, n_comp_diesel, C_comp_diesel, n_exp_diesel, C_exp_diesel, 'Diesel', ...
                                     P_petrol_1400, V_petrol_1400, states_petrol, n_comp_petrol, C_comp_petrol, n_exp_petrol, C_exp_petrol, 'Petrol', ...
                                     1400, 'figures/dual_engine_1400RPM', save_plots);
end

%% --- Save Results to Excel ---

filename = 'results_tables.xlsx';
writetable(engine_data_table, filename, 'Sheet', 'Raw Data');
writetable(engine_performance_table, filename, 'Sheet', 'Performance');
writetable(heat_balance_table, filename, 'Sheet', 'Heat Balance');
writetable(PV_analysis_table_combined, filename, 'Sheet', 'PV Analysis');
writetable(PV_states_table, filename, 'Sheet', 'PV States');

%% Function Definitions

function [performance_row, heat_balance_row] = calculate_data_row(row, V_d, density_fuel, density_air, Q_HV, c_p_water, c_p_exhaust, n_r_strokes)
    % Extract data
    nom_rpm = row.NominalSpeed;
    engine_load = row.EngineLoad;
    rpm = row.MeasuredSpeed;
    torque = row.Torque;
    power = row.Power;
    elevation_diff = row.ElevationDiff;
    vol_flow_fuel_lpm = row.FuelRate;
    vol_flow_water_lpm = row.CoolRate;
    temp_water_in = row.WaterIn;
    temp_water_out = row.WaterOut;
    temp_cal_out = row.CalOut;
    temp_exh_in = row.ExhIn;

    % Calculate Engine Performance
    [performance_row, mass_flow_fuel, mass_flow_air] = ...
        calculate_engine_performance(nom_rpm, engine_load, rpm, torque, power, ...
                                     elevation_diff, vol_flow_fuel_lpm, ...
                                     density_fuel, density_air, Q_HV, V_d, n_r_strokes);

    % Convert flow rates to kg/s
    mass_flow_water = (vol_flow_water_lpm / 60000) * 1000; % (L/min -> m^3/s) * (kg/m^3)

    % Calculate Heat Balance
    heat_balance_row = calculate_heat_balance(nom_rpm, engine_load, mass_flow_air, mass_flow_fuel, mass_flow_water, c_p_water, c_p_exhaust, temp_water_in, temp_water_out, temp_exh_in, temp_cal_out, Q_HV, power);
end

function [row, mass_flow_fuel, mass_flow_air] = calculate_engine_performance(nom_rpm, engine_load, rpm, torque, power, elevation_diff, vol_flow_fuel_lpm, density_fuel, density_air, Q_HV, V_d, n_r_strokes)
    
    vol_flow_fuel_m3s = (vol_flow_fuel_lpm / 60000); % m^3 s^-1
    mass_flow_fuel = density_fuel * vol_flow_fuel_m3s; % kg s^-1
    
    % Air Flow
    [mass_flow_air, vol_flow_air] = mass_and_volume_flow_rate_air(elevation_diff, density_air);
    
    % Performance Metrics
    BMEP = (2 * pi * n_r_strokes * torque / V_d) / 1e6; % MPa
    SFC = (mass_flow_fuel / power) * 3600 * 1000; % g/kW-hr
    air_fuel_ratio = mass_flow_air / mass_flow_fuel;
    efficiency_vol = (mass_flow_air * n_r_strokes) / (density_air * V_d * (rpm / 60)) * 100;
    efficiency_thermal = power / (mass_flow_fuel * Q_HV) * 100;

    % Create the result row
    row = [nom_rpm, engine_load, rpm, power, torque, BMEP, SFC, air_fuel_ratio, ...
           efficiency_vol, efficiency_thermal, vol_flow_fuel_m3s, mass_flow_fuel, ...
           vol_flow_air, mass_flow_air];
end

function row = calculate_heat_balance(nom_rpm, engine_load, mass_flow_air, mass_flow_fuel, mass_flow_water, c_p_water, c_p_exhaust, T_water_in, T_water_out, T_exh_in, T_cal_out, Q_HV, brake_power)
    % Heat balance calculations. The calculated values are in power.
    Q_fuel = mass_flow_fuel * Q_HV; % kW
    Q_cooling = mass_flow_water * c_p_water * (T_water_out - T_water_in); % kW
    mass_flow_exhaust = mass_flow_air + mass_flow_fuel; % kg/s
    Q_exhaust = mass_flow_exhaust * c_p_exhaust * (T_exh_in - T_cal_out); % kW
    Q_other = Q_fuel - Q_cooling - Q_exhaust - brake_power; % kW (unaccounted heat)

    row = [nom_rpm, engine_load, Q_fuel, Q_cooling, Q_exhaust, brake_power, Q_other, mass_flow_exhaust, mass_flow_water];
end

function [mass_flow_rate, vol_flow_rate] = ...
    mass_and_volume_flow_rate_air(elevation_diff, density_air)
    % Calculates the mass [kg s^-1] and volume [m^3 s^-1] flow rate into the 
    % engine using an orifice meter.
    % elevation_diff: Difference in height of the manometer fluid [mm]
    % density_air: Density of the air outside the manometer [kg m^-3]
    
    C_0 = 0.6;          % Correction factor for air flowing through orifice
    orifice_diameter = 0.03; % Orifice diameter [m]
    density_man = 1000; % Density of the manometer fluid (water) [kg m^-3]
    g = 9.81;           % Gravitational acceleration [m s^-2]
    
    A_0 = pi * orifice_diameter^2 / 4; % Orifice Area [m^2]
    
    % Pressure difference [Pa]
    pressure_diff = (density_man - density_air) * g * (elevation_diff / 1000); 
    
    % Volumetric Flow Rate [m^3/s]
    vol_flow_rate = C_0 * A_0 * sqrt(2 * pressure_diff / density_air);
    
    % Mass Flow Rate [kg/s]
    mass_flow_rate = vol_flow_rate * density_air;
end

function plot_engine_performance(result_table, filename, save_plots)
    % Function to plot engine performance parameters vs RPM for two loads across three RPMs.
    %
    % Assumes the following data structure in result_table (rows 1-6):
    % Rows 1, 3, 5 are for Load A (e.g., 50% load) at three different RPMs.
    % Rows 2, 4, 6 are for Load B (e.g., 100% load) at three different RPMs.
    
    % --- Data Extraction for Multi-RPM Plot ---
    
    % Data for Load A (e.g., 50% load): Rows 1, 3, 5
    rows_A = [1, 3, 5];
    data_A = result_table(rows_A, :);
    
    % Data for Load B (e.g., 100% load): Rows 2, 4, 6
    rows_B = [2, 4, 6];
    data_B = result_table(rows_B, :);
    
    % Extract vectors for Load A
    rpm_A = data_A.Measured_RPM;
    torque_A = data_A.Torque_Nm;
    power_A = data_A.Power_kW;
    BMEP_A = data_A.BMEP_MPa;
    SPC_A = data_A.SFC_g_kW_hr;
    air_fuel_A = data_A.Air_Fuel_Ratio;
    efficiency_vol_A = data_A.Vol_Efficiency_pct;
    efficiency_thermal_A = data_A.Thermal_Efficiency_pct;
    
    % Extract vectors for Load B
    rpm_B = data_B.Measured_RPM;
    torque_B = data_B.Torque_Nm;
    power_B = data_B.Power_kW;
    BMEP_B = data_B.BMEP_MPa;
    SPC_B = data_B.SFC_g_kW_hr;
    air_fuel_B = data_B.Air_Fuel_Ratio;
    efficiency_vol_B = data_B.Vol_Efficiency_pct;
    efficiency_thermal_B = data_B.Thermal_Efficiency_pct;
    
    % --- Plotting Setup ---
    
    % Set font properties (Reused from original)
    font_name = 'Times New Roman';
    font_size_label = 14; 
    font_size_title = 16;
    font_size_legend = 20;
    line_thickness = 1.0; % Reusing marker_thickness as line thickness
    marker_size_plot = 8; % Appropriate size for 'plot' markers
    
    % Create subplots
    figure('Units', 'Inches', 'Position', [1, 1, 10, 12]);
    
    % Title (Changed to general title)
    title_string = 'Engine Performance vs RPM at Varying Loads';
    sgtitle(title_string, 'FontName', font_name, 'FontSize', font_size_title, 'FontWeight', 'bold');
    
    % --- Subplots: Changed to 'plot' for line curves ---
    
    % 1. Torque vs RPM
    subplot(4, 2, 1);
    % Changed from scatter to plot with line and markers ('-b*' for Load A)
    plot(rpm_A, torque_A, '-b*', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    hold on;
    % Changed from scatter to plot with line and markers ('-ro' for Load B)
    plot(rpm_B, torque_B, '-ro', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    title('Torque vs RPM', 'FontName', font_name, 'FontSize', font_size_label);
    xlabel('RPM', 'FontName', font_name, 'FontSize', font_size_label);
    ylabel('Torque (Nm)', 'FontName', font_name, 'FontSize', font_size_label);
    grid on;
    set(gca, 'Box', 'on');
    hold off;
    
    % 2. Brake Power vs RPM
    subplot(4, 2, 2);
    plot(rpm_A, power_A, '-b*', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    hold on;
    plot(rpm_B, power_B, '-ro', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    title('Brake Power vs RPM', 'FontName', font_name, 'FontSize', font_size_label);
    xlabel('RPM', 'FontName', font_name, 'FontSize', font_size_label);
    ylabel('Power (kW)', 'FontName', font_name, 'FontSize', font_size_label);
    grid on;
    set(gca, 'Box', 'on');
    hold off;
    
    % 3. BMEP vs RPM
    subplot(4, 2, 3);
    plot(rpm_A, BMEP_A, '-b*', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    hold on;
    plot(rpm_B, BMEP_B, '-ro', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    title('BMEP vs RPM', 'FontName', font_name, 'FontSize', font_size_label);
    xlabel('RPM', 'FontName', font_name, 'FontSize', font_size_label);
    ylabel('BMEP (MPa)', 'FontName', font_name, 'FontSize', font_size_label);
    grid on;
    set(gca, 'Box', 'on');
    hold off;
    
    % 4. Specific Fuel Consumption vs RPM
    subplot(4, 2, 4);
    plot(rpm_A, SPC_A, '-b*', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    hold on;
    plot(rpm_B, SPC_B, '-ro', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    title('Specific Fuel Consumption vs RPM', 'FontName', font_name, 'FontSize', font_size_label);
    xlabel('RPM', 'FontName', font_name, 'FontSize', font_size_label);
    ylabel('SFC (g/kW-hr)', 'FontName', font_name, 'FontSize', font_size_label);
    grid on;
    set(gca, 'Box', 'on');
    hold off;
    
    % 5. Air-Fuel Ratio vs RPM
    subplot(4, 2, 5);
    plot(rpm_A, air_fuel_A, '-b*', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    hold on;
    plot(rpm_B, air_fuel_B, '-ro', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    title('Air-Fuel Ratio vs RPM', 'FontName', font_name, 'FontSize', font_size_label);
    xlabel('RPM', 'FontName', font_name, 'FontSize', font_size_label);
    ylabel('Air-Fuel Ratio', 'FontName', font_name, 'FontSize', font_size_label);
    grid on;
    set(gca, 'Box', 'on');
    hold off;
    
    % 6. Volumetric Efficiency vs RPM
    subplot(4, 2, 6);
    plot(rpm_A, efficiency_vol_A, '-b*', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    hold on;
    plot(rpm_B, efficiency_vol_B, '-ro', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    title('Volumetric Efficiency vs RPM', 'FontName', font_name, 'FontSize', font_size_label);
    xlabel('RPM', 'FontName', font_name, 'FontSize', font_size_label);
    ylabel('Volumetric Efficiency (%)', 'FontName', font_name, 'FontSize', font_size_label);
    ylim([0 100]); 
    grid on;
    set(gca, 'Box', 'on');
    hold off;
    
    % 7. Thermal Efficiency vs RPM
    subplot(4, 2, 7);
    h1 = plot(rpm_A, efficiency_thermal_A, '-b*', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    hold on;
    h2 = plot(rpm_B, efficiency_thermal_B, '-ro', 'LineWidth', line_thickness, 'MarkerSize', marker_size_plot);
    title('Thermal Efficiency vs RPM', 'FontName', font_name, 'FontSize', font_size_label);
    xlabel('RPM', 'FontName', font_name, 'FontSize', font_size_label);
    ylabel('Thermal Efficiency (%)', 'FontName', font_name, 'FontSize', font_size_label);
    ylim([0 100]); 
    grid on;
    set(gca, 'Box', 'on');
    hold off;
    
    % 8. Legend Plot
    legend_plot = subplot(4, 2, 8);
    position_legend_plot = get(legend_plot,'position');
    axis(legend_plot, 'off'); % Turn off axes for the legend space
    
    % Call the legend using the handles (h1, h2) from subplot 7.
    legend_labels = {sprintf('%i%% Load', 50), sprintf('%i%% Load', 100)};
    legend_all_plots = legend(legend_plot, [h1, h2], legend_labels, ...
        'FontName', font_name, 'FontSize', font_size_legend, 'Box', 'on', 'EdgeColor', 'k', 'Location', 'NorthWest');
    
    % Adjust the position of the legend to fill the empty subplot space
    set(legend_all_plots, 'position', position_legend_plot); 
    set(get(legend_all_plots, 'Title'), 'String', 'Load');

    % Save the figure
    if save_plots
        figure_name = sprintf('%s.tif', filename);
        exportgraphics(gcf, figure_name, "Resolution", 600)
    end
end

function plot_heat_balance_pie_chart(row_index, result_table, filename, save_plots)
    % Set font properties for pie chart labels
    pie_label_font_name = 'Times New Roman';
    pie_label_font_size = 14; 
    pie_title_font_size = 16;
    
    % Extract the relevant data row
    row = result_table(row_index, :);
    
    % Energy output components (Losses and Useful Work)
    output_data = [row.Q_cooling_kW, row.Q_exhaust_kW, row.Brake_Power_kW, row.Q_other_kW];
    
    % Calculate percentages for labels
    fuel_total = row.Q_fuel_kW;
    Q_cool_pct = (row.Q_cooling_kW / fuel_total) * 100;
    Q_exh_pct = (row.Q_exhaust_kW / fuel_total) * 100;
    Brake_pct = (row.Brake_Power_kW / fuel_total) * 100;
    Q_other_pct = (row.Q_other_kW / fuel_total) * 100;
    
    output_labels = {
        sprintf('\\bf{Q_{Cooling}}\n%.1f%%\n%.1f kW', Q_cool_pct, row.Q_cooling_kW), ...
        sprintf('\\bf{Q_{Exhaust}}\n%.1f%%\n%.1f kW', Q_exh_pct, row.Q_exhaust_kW), ...
        sprintf('\\bf{W_{Brake}}\n%.1f%%\n%.1f kW', Brake_pct, row.Brake_Power_kW), ...
        sprintf('\\bf{Q_{Other}}\n%.1f%%\n%.1f kW', Q_other_pct, row.Q_other_kW)}; 
    
    figure;
    h = pie(output_data, output_labels);
    
    % Set font properties for pie chart labels
    for k = 1:length(h)
        if isgraphics(h(k), 'text')
            h(k).FontName = pie_label_font_name;
            h(k).FontSize = pie_label_font_size;
        end
    end
    
    % Title 
    title_str = sprintf('Heat Balance at %i RPM and %i%% Engine Load\nQ_{Fuel} = %.2f kW', ...
        row.Nominal_RPM, row.Engine_Load_pct, row.Q_fuel_kW);
    title(title_str, ...
        'FontSize', pie_title_font_size, ...
        'FontWeight', 'bold', ...
        'FontName', pie_label_font_name);
    
    axis equal; % Ensure the pie chart is round
    
    % Set figure size
    set(gcf, 'Units', 'Inches', 'Position', [1, 1, 10, 10]);

    if save_plots
        % Save the figure
        figure_name = sprintf('%s.tif', filename);
        exportgraphics(gca, figure_name, "Resolution", 300)
    end
end

function [n_comp, C_comp, n_exp, C_exp, states] = polytropic_coefficient_fitting(P_data, V_data, fitted_points, nominal_rpm, T1, k_air, draw_plots, filename, save_plots)
    % Plot parameters
    markerSize = 100; % Size of marker on plot

    % Visual Line Fitting for Polytropic Coefficients
    % Adjust these variables until the lines fit the data.
    % Use the cursor on the plot to find coordinates for volume and pressure of
    % the start and end of the linear sections of data.
    
    P1_comp = fitted_points(1); V1_comp = fitted_points(2);
    P2_comp = fitted_points(3); V2_comp = fitted_points(4);
    P1_exp = fitted_points(5); V1_exp = fitted_points(6);
    P2_exp = fitted_points(7); V2_exp = fitted_points(8);
    
    % Calculations for polytropic coefficients
    % The log-log plot has linear sections for the expansion and compression,
    % so by fitting a linear curve visually by setting points the slope of that
    % curve can be calculated, giving the coefficients for both.
    
    % Calculation for Compression Line
    % This calculates the slope of the line (n, ie the polytropic coefficient)
    n_comp = -1 * (log10(P2_comp) - log10(P1_comp)) / (log10(V2_comp) - log10(V1_comp));
    % Calculate the constant
    C_comp = P1_comp * V1_comp^n_comp;
    
    % Calculation for Expansion Line
    n_exp = -1 * (log10(P2_exp) - log10(P1_exp)) / (log10(V2_exp) - log10(V1_exp));
    C_exp = P1_exp * V1_exp^n_exp; 
    
    % Display the polytropic coefficients
    %{
    fprintf('\nPolytropic Coefficients: %i RPM: %s\n', nominal_rpm, filename);
    fprintf('n-Compression = %.3f, C-Compression = %.3f\n', n_comp, C_comp);
    fprintf('n-Expansion   = %.3f, C-Expansion   = %.3f\n', n_exp, C_exp);
    %}

    % Get the points for each state
    V1 = V1_comp;
    P1 = C_comp / V1^n_comp;
    
    % Get the indices for Pressure values greater than half the max pressure,
    % to filter out the lower exhaust cycle.
    V_min_indicies = find(P_data > (max(P_data)/2));
    % Get the volume values associated with the pressure indicies
    V_min_candidates = V_data(V_min_indicies);
    % Find the minimum volume and its index
    V2 = min(V_min_candidates);
    P2 = C_comp / V2^n_comp;
    
    P3 = P2;
    V3 = (C_exp / P3)^(1 / n_exp);
    
    V4 = V1;
    P4 = C_exp / V4^n_exp;

    compression_ratio = V1 / V2;
    T2 = T1 * compression_ratio^(k_air - 1);
    T3 = T2 * (V3 / V2);
    T4 = T3 * (1/ compression_ratio)^(k_air - 1);

    states = [P1, V1, T1, P2, V2, T2, P3, V3, T3, P4, V4, T4];

    % -- Plotting --
    if draw_plots
        % Define font size variables for title, legend, and axis labels
        title_font_size = 16;
        legend_font_size = 12;
        axis_label_font_size = 14;

        % Define font name variable
        font_name = 'Times New Roman';

        % Calculate x and y values. P = C * V^n
        V_fit_comp_log = logspace(log10(min(V_data)), log10(max(V_data)), 100); 
        P_fit_comp_log = C_comp * V_fit_comp_log.^-n_comp;
        V_fit_exp_log = logspace(log10(min(V_data)), log10(max(V_data)), 100); 
        P_fit_exp_log = C_exp * V_fit_exp_log.^-n_exp;
      
        figure;
        % Plot experimental data
        scatter(V_data, P_data, 'filled', 'DisplayName', 'P-V Data');
        hold on;
        
        % Plot fitted compression line
        plot(V_fit_comp_log, P_fit_comp_log, 'r-', 'LineWidth', 2, ...
            'DisplayName', ['Fit 1 (n = ', num2str(n_comp, '%.3f'), ')']);
        % Scatter plot for the compression stroke points
        scatter([V1_comp, V2_comp], [P1_comp, P2_comp], markerSize, 'ro', 'filled', 'DisplayName', 'Compression Stroke Points');
        
        % Plot fitted expansion line
        plot(V_fit_exp_log, P_fit_exp_log, 'g-', 'LineWidth', 2, ...
            'DisplayName', ['Fit 2 (n = ', num2str(n_exp, '%.3f'), ')']);
        % Scatter plot for the expansion stroke points
        scatter([V1_exp, V2_exp], [P1_exp, P2_exp], markerSize, 'go', 'filled', 'DisplayName', 'Expansion Stroke Points');
        
        title_string = sprintf('Log-Log P-V Diagram at %i RPM and 100%% Load', nominal_rpm);
        title(title_string, 'FontSize', title_font_size, 'FontName', font_name);
        xlabel('Volume (m^3)', 'FontSize', axis_label_font_size, 'FontName', font_name);
        ylabel('Pressure (kPa)', 'FontSize', axis_label_font_size, 'FontName', font_name);
        grid on;
        set(gca, 'XScale', 'log', 'YScale', 'log'); % Set log-log scale
        
        legend('Location','northeast', 'FontSize', legend_font_size, 'FontName', font_name);
        ylim([min(P_data)*0.9 max(P_data)*1.1]);
        xlim([min(V_data)*0.9 max(V_data)*1.1]);
        hold off;
    
        % Set the size of the figure
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, 10, 6]);
    
        if save_plots
            % Save the figure
            figure_name = sprintf('%s.tif', filename);
            exportgraphics(gca, figure_name, "Resolution", 300)
        end
    end
end

function plot_fitted_PV_curve(P_data, V_data, n_comp, C_comp, n_exp, C_exp, states, nominal_rpm, filename, save_plots)
    P1 = states(1);  V1 = states(2);  T1 = states(3);
    P2 = states(4);  V2 = states(5);  T2 = states(6);
    P3 = states(7);  V3 = states(8);  T3 = states(9);
    P4 = states(10); V4 = states(11); T4 = states(12);
    
    marker_size_PV_data = 10; % Size of marker for states
    marker_size_states = 50; % Size of marker for states
    
    fitted_curve_line_width = 2;
    fitted_curve_style = 'g--';
    
    title_font_size = 16;
    label_font_size = 14;
    legend_font_size = 12;
    
    figure;
    scatter(V_data, P_data, marker_size_PV_data, 'k', 'filled');
    hold on;
    
    % State 1
    scatter(V1, P1, marker_size_states, 'r', 'filled');
    % State 2
    scatter(V2, P2, marker_size_states, 'g', 'filled');
    % State 3
    scatter(V3, P3, marker_size_states, 'm', 'filled');
    % State 4
    scatter(V4, P4, marker_size_states, 'cyan', 'filled');
    
    % Fitted PV Curve for Compression and Expansion
    V_fit_comp = linspace(V2, V1, 100);
    P_fit_comp = C_comp * V_fit_comp.^-n_comp;
    V_fit_exp = linspace(V3, V4, 100);
    P_fit_exp = C_exp * V_fit_exp.^-n_exp;
    
    % Combine the compression and expansion curves into a single plot
    plot(V_fit_comp, P_fit_comp, fitted_curve_style, 'LineWidth', fitted_curve_line_width);
    hold on;
    plot(V_fit_exp, P_fit_exp, fitted_curve_style, 'LineWidth', fitted_curve_line_width); 
    
    % Connecting points 4 to 1
    plot([V4, V1], [P4, P1], fitted_curve_style, 'LineWidth', fitted_curve_line_width);
    
    % Connecting points 2 to 3
    plot([V2, V3], [P2, P3], fitted_curve_style, 'LineWidth', fitted_curve_line_width);
    
    % Update legend to include all markers and curves
    legend('P-V Data', 'State 1', 'State 2', 'State 3', 'State 4', 'Fitted Diesel Cycle', 'Location', 'northeast', 'FontSize', legend_font_size);
    xlabel('Volume (m^3)', 'FontName', 'Times New Roman', 'FontSize', label_font_size);
    ylabel('Pressure (kPa)', 'FontName', 'Times New Roman', 'FontSize', label_font_size);
    title_string = sprintf('P-V Diagram with Fitted Diesel Cycle at %i RPM', nominal_rpm);
    title(title_string, 'FontSize', title_font_size);
    set(gca, 'FontName', 'Times New Roman'); % Set font for axes
    grid on;
    hold off;
    
    % Set the size of the figure
    set(gcf, 'Units', 'Inches', 'Position', [1, 1, 10, 6]);

    if save_plots
        % Save the figure
        figure_name = sprintf('%s.tif', filename);
        exportgraphics(gca, figure_name, "Resolution", 300)
    end
end

function row = PV_analysis(V, P, rpm, power, n_r_strokes, number_cylinders, V_d, mass_flow_air, mass_flow_fuel, Q_HV, c_p_air, k_air, states)
    P1 = states(1);  V1 = states(2);  T1 = states(3);
    P2 = states(4);  V2 = states(5);  T2 = states(6);
    P3 = states(7);  V3 = states(8);  T3 = states(9);
    P4 = states(10); V4 = states(11); T4 = states(12);
    
    compression_ratio = V1 / V2;
    cut_off_ratio = V3 / V2;
    T2 = T1 * compression_ratio^(k_air - 1);
    T3 = T2 * (V3 / V2);
    
    indicated_work = trapz(V, P); % kJ

    indicated_power = (indicated_work * (rpm / 60) * number_cylinders) / n_r_strokes; % kW

    mechanical_efficiency = (power / indicated_power) * 100;
    
    IMEP = (indicated_work * number_cylinders / V_d) / 1000; %MPa
    
    mass_fuel = (mass_flow_fuel * n_r_strokes) / ((rpm / 60) * number_cylinders); % kg
    mass_air = (mass_flow_air * n_r_strokes) / ((rpm / 60) * number_cylinders); % kg
    energy_fuel = mass_fuel * Q_HV; % kJ
    energy_in = mass_air * c_p_air * (T3 - T2); % kJ
    efficiency = (energy_in / energy_fuel) * 100;

    row = [indicated_work, indicated_power, mechanical_efficiency, IMEP, efficiency, cut_off_ratio, energy_in, energy_fuel, mass_fuel, mass_air];
end

function [PV_results_row, states_row] = process_pv_data(V_data, P_data, fitted_points, perf_table, row_index, n_r_strokes, number_cylinders, V_d, Q_HV, c_p_air, k_air, T1, draw_plots, save_plots)
    
    % Extract data from Engine Performance Table
    perf_row = perf_table(row_index, :);
    nominal_rpm = perf_row.Nominal_RPM;
    measured_rpm = perf_row.Measured_RPM;
    brake_power = perf_row.Power_kW;
    mass_flow_air = perf_row.Mass_Flow_Air_kg_s;
    mass_flow_fuel = perf_row.Mass_Flow_Fuel_kg_s;
    
    % Polytropic Coefficient Fitting
    [n_comp, C_comp, n_exp, C_exp, states] = polytropic_coefficient_fitting(...
        P_data, V_data, fitted_points, nominal_rpm, T1, k_air, draw_plots,...
        sprintf('figures/polytropic_fitting_curve_%i', nominal_rpm), save_plots);
    
    % Plot Fitted PV Curve
    if draw_plots
        plot_fitted_PV_curve(P_data, V_data, n_comp, C_comp, n_exp, C_exp, states, nominal_rpm, sprintf('figures/diesel_cycle_%i', nominal_rpm), save_plots);
    end
    
    % PV Analysis Calculation
    PV_analysis_results = PV_analysis(V_data, P_data, measured_rpm, brake_power, ...
        n_r_strokes, number_cylinders, V_d, mass_flow_air, mass_flow_fuel, Q_HV, ...
        c_p_air, k_air, states);
    
    PV_results_row = [nominal_rpm, measured_rpm, PV_analysis_results];
    states_row = [nominal_rpm, n_comp, C_comp, n_exp, C_exp, states];
end

function plot_fitted_PV_curve_dual_engine(P_data1, V_data1, states1, n_comp1, C_comp1, n_exp1, C_exp1, engine1_type, P_data2, V_data2, states2, n_comp2, C_comp2, n_exp2, C_exp2, engine2_type, nominal_rpm, filename, save_plots)
% Plots two independent sets of P-V data and their corresponding fitted cycles.
%
% Inputs:
% P_data1, V_data1, states1, n_comp1, C_comp1, n_exp1, C_exp1: Parameters for Engine 1.
% P_data2, V_data2, states2, n_comp2, C_comp2, n_exp2, C_exp2: Parameters for Engine 2.
% nominal_rpm: Engine speed for the title.
% filename: Base name for saving the plot file.
% save_plots: Boolean to indicate whether to save the plot.

    % --- Engine 1 Parameters (Blue/Black) ---
    engine1_color_data = [0 0 0];       % Black data points
    engine1_color_fit  = 'g--';         % Green fitted curve
    engine1_states = states1;
    
    % --- Engine 2 Parameters (Red) ---
    engine2_color_data = [0.5 0.5 0.5]; % Gray data points
    engine2_color_fit  = 'r--';         % Red fitted curve
    engine2_states = states2;

    % --- Plot Style Definitions (Common) ---
    marker_size_PV_data = 10;
    fitted_curve_line_width = 2;
    
    title_font_size = 16;
    label_font_size = 14;
    legend_font_size = 12;
    
    figure;
    
    % --- Plot Engine 1 ---

    % Unpack State Points for Engine 1
    P1_1 = engine1_states(1);  V1_1 = engine1_states(2);  T1_1 = engine1_states(3);
    P2_1 = engine1_states(4);  V2_1 = engine1_states(5);  T2_1 = engine1_states(6);
    P3_1 = engine1_states(7);  V3_1 = engine1_states(8);  T3_1 = engine1_states(9);
    P4_1 = engine1_states(10); V4_1 = engine1_states(11); T4_1 = engine1_states(12);
    % Plot Data Set 1
    engine1_PV_legend_string = sprintf('%s P-V Data', engine1_type);
    scatter(V_data1, P_data1, marker_size_PV_data, engine1_color_data, 'filled', 'DisplayName', engine1_PV_legend_string);
    hold on;

    % Fitted PV Curve for Engine 1
    V_fit_comp1 = linspace(V2_1, V1_1, 100);
    P_fit_comp1 = C_comp1 * V_fit_comp1.^-n_comp1;
    V_fit_exp1  = linspace(V3_1, V4_1, 100);
    P_fit_exp1  = C_exp1 * V_fit_exp1.^-n_exp1;
    
    % Plot Fitted Diesel Cycle 1 (Use a solid line for the primary plot)
    engine1_fitted_legend_string = sprintf('%s Fitted Cycle', engine1_type);
    plot(V_fit_comp1, P_fit_comp1, engine1_color_fit, 'LineWidth', fitted_curve_line_width, 'DisplayName', engine1_fitted_legend_string);
    plot(V_fit_exp1, P_fit_exp1, engine1_color_fit, 'LineWidth', fitted_curve_line_width, 'HandleVisibility', 'off'); 
    
    % Connecting points for Engine 1 (Use 'HandleVisibility', 'off' to keep off legend)
    plot([V4_1, V1_1], [P4_1, P1_1], engine1_color_fit, 'LineWidth', fitted_curve_line_width, 'HandleVisibility', 'off');
    plot([V2_1, V3_1], [P2_1, P3_1], engine1_color_fit, 'LineWidth', fitted_curve_line_width, 'HandleVisibility', 'off');
    
    % --- PLot Engine 2 ---
    
    % Unpack State Points for Engine 2
    P1_2 = engine2_states(1);  V1_2 = engine2_states(2);
    P2_2 = engine2_states(4);  V2_2 = engine2_states(5);
    P3_2 = engine2_states(7);  V3_2 = V2_2; % Cheating since its a petrol cycle, also manually set P3 on the plot
    P4_2 = engine2_states(10); V4_2 = engine2_states(11);
    % Plot Data Set 2
    engine2_PV_legend_string = sprintf('%s P-V Data', engine2_type);
    scatter(V_data2, P_data2, marker_size_PV_data, engine2_color_data, 's', 'DisplayName', engine2_PV_legend_string); % Use 's' for square markers
    
    % Fitted PV Curve for Engine 2
    V_fit_comp2 = linspace(V2_2, V1_2, 100);
    P_fit_comp2 = C_comp2 * V_fit_comp2.^-n_comp2;
    V_fit_exp2  = linspace(V3_2, V4_2, 100);
    P_fit_exp2  = C_exp2 * V_fit_exp2.^-n_exp2;
    
    % Plot Fitted Diesel Cycle 2 (Use a dashed line for distinction)
    engine2_fitted_legend_string = sprintf('%s Fitted Cycle', engine2_type);
    plot(V_fit_comp2, P_fit_comp2, engine2_color_fit, 'LineWidth', fitted_curve_line_width, 'DisplayName', engine2_fitted_legend_string);
    plot(V_fit_exp2, P_fit_exp2, engine2_color_fit, 'LineWidth', fitted_curve_line_width, 'HandleVisibility', 'off'); 
    
    % Connecting points for Engine 2
    plot([V4_2, V1_2], [P4_2, P1_2], engine2_color_fit, 'LineWidth', fitted_curve_line_width, 'HandleVisibility', 'off');
    plot([V2_2, V3_2], [P2_2, 5250], engine2_color_fit, 'LineWidth', fitted_curve_line_width, 'HandleVisibility', 'off');

    % --- Finalise Plot ---

    % Update legend: It automatically picks up items with 'DisplayName'
    legend('Location', 'northeast', 'FontSize', legend_font_size);

    xlabel('Volume (m^3)', 'FontName', 'Times New Roman', 'FontSize', label_font_size);
    ylabel('Pressure (kPa)', 'FontName', 'Times New Roman', 'FontSize', label_font_size);
    title_string = sprintf('Diesel vs. Petrol Engine P-V Comparison at %i RPM', nominal_rpm);
    title(title_string, 'FontSize', title_font_size);
    set(gca, 'FontName', 'Times New Roman');
    grid on;
    hold off;
    
    % Set the size of the figure
    set(gcf, 'Units', 'Inches', 'Position', [1, 1, 10, 6]);
    
    % Save Plot
    if save_plots
        figure_name = sprintf('%s.tif', filename);
        exportgraphics(gca, figure_name, "Resolution", 300)
    end
end

