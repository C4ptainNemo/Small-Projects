%% EGB 323 Practical
% This script will calculate and output a table for the EGB323 Practical data.
% Written by Adam Telfer.
%
% The fluid properties section applies to all 4 pipe sections.
% Individual pipe sections have specific data in them to change.
% For the arrays, ensure they are the same length for each pipe (nominally
%   5 data points).
% Units are in SI unless specified, eg. mmHg and L/min which come from
%   practical measurements.

%---------- Equations ----------%
% This is the key equations applied throughout, friction factors are contained
%   within their functions.
%
% Pa = 133.322 * mmHg
%
% P = density_water * g * head_loss => headl_loss = P / (density_water * g)
%
% Water Velocity:   flow_rate = (flow_rate_Lmin / 60) / 1000; % Convert L/min to m^3/s
%                   area = (pi / 4) * d^2;
%                   velocity = flow_rate / area;
%
% Reynold Number: Re = (density_water * velocity * pipe_diameter) / dynamic_viscosity
%
% Relative Pipe Roughness (epsilon/D): rel_roughness = pipe_roughness / pipe_diameter;
%
% Percentage Error = abs(experimental - theoretical) / theoretical * 100;

%---------- Friction Factors ----------%
% Three different friction factors have been used to compute the
%   theoretical head loss within the pipe.
% Moody: Values taken by calculating the relative roughness of the pipe and
%   Reynolds number of the flow and looking it up on a Moody chart
% Colebrook-White: MATLAB can solve a non-linear function, check function
%   description for more details.
% Haaland: The Haaland friction factor equation was used (again check
%   function for implementation). It is only valid for Re > 4000 (though
%   this isn't an issue for most real flows).

% ----- Clear ----- %
clear;      % Clear Workspace
close all;  % Close all figures
clc;        % Clear Command Window

%---------- Fluid Properties ----------%
% Set the two temp values, run script and use printed out value for density and viscosity
% from https://webbook.nist.gov/chemistry/fluid/
% Pressure was assumed at 1 atm
water_temp_start = 25.8; % C, water temperature at start of experiment
water_temp_end = 26.6; % C, water temperature at end of experiment
water_temp_avg = (water_temp_start + water_temp_end) / 2;

water_density = 996.73; % kg/m^3 at average temp (26.2)
dynamic_viscosity_water = 866.21e-6; % Pa*s at average temp (26.2)

fprintf(1, 'T_start = %g C, T_end = %g C, T_avg = %g C\n', ...
            water_temp_start, water_temp_end, water_temp_avg)
fprintf(1, ['Water Density = %g kg/m^3, ' ...
            'Dynamic Viscosity = %.3e Pa*s\n' ...
            'Values derived from "https://webbook.nist.gov/chemistry/fluid/" ' ...
            'using average temperature and 1 atm\n'], ...
            water_density, dynamic_viscosity_water)

%---------- Pipe 1 ----------%
fprintf(1, '\nPipe 1\n')
pipe_diameter = 10.92e-3; % 10.92 mm
pipe_length = 1.22; % m
pipe_roughness = 1.5e-6; % 1.5e-3 mm
flow_rate = [8, 16, 24, 32, 34]; % L/min
manometer_reading_mmHg = [22.5, 67.8, 134.6, 232, 255]; % mmHg
% Run script, get Re, change default values
moody_friction_factor = [0.028, 0.023, 0.021, 0.020, 0.019];

pipe1_results = pipe_results(flow_rate, ...
                                       manometer_reading_mmHg, ...
                                       dynamic_viscosity_water, ...
                                       pipe_diameter, ...
                                       pipe_roughness, ...
                                       pipe_length, ...
                                       water_density, ...
                                       moody_friction_factor);

disp(struct2table(pipe1_results));
plot_flow_vs_headloss(pipe1_results, 'Pipe 1', 4);
writetable(struct2table(pipe1_results), 'pipe1_results.csv');

%---------- Pipe 2 ----------%
fprintf(1, '\nPipe 2\n')
pipe_diameter = 13.84e-3; % 13.84 mm
pipe_length = 1.22; % m
pipe_roughness = 1.5e-6; % 1.5e-3 mm
flow_rate = [8, 16, 24, 32, 40]; % L/min
manometer_reading_mmHg = [8.5, 24.0, 48.4, 80.0, 120.0]; % mmHg
% Run script, get Re, change default values
moody_friction_factor = [0.030, 0.025, 0.022, 0.021, 0.02]; 

pipe2_results = pipe_results(flow_rate, ...
                                       manometer_reading_mmHg, ...
                                       dynamic_viscosity_water, ...
                                       pipe_diameter, ...
                                       pipe_roughness, ...
                                       pipe_length, ...
                                       water_density, ...
                                       moody_friction_factor);

disp(struct2table(pipe2_results));
plot_flow_vs_headloss(pipe2_results, 'Pipe 2', 4);
writetable(struct2table(pipe2_results), 'pipe2_results.csv');

%---------- Pipe 3 ----------%
fprintf(1, '\nPipe 3\n')
pipe_diameter = 19.94e-3; % 19.94 mm
pipe_length = 1.22; % m
pipe_roughness = 1.5e-6; % 1.5e-3 mm
flow_rate = [8, 16, 24, 32, 40]; % L/min
manometer_reading_mmHg = [0, 0, 7.8, 15.3, 22.2]; % mmHg
% Run script, get Re, change default values
moody_friction_factor = [0.034, 0.028, 0.025, 0.023, 0.021]; 

pipe3_results = pipe_results(flow_rate, ...
                                       manometer_reading_mmHg, ...
                                       dynamic_viscosity_water, ...
                                       pipe_diameter, ...
                                       pipe_roughness, ...
                                       pipe_length, ...
                                       water_density, ...
                                       moody_friction_factor);

disp(struct2table(pipe3_results));
plot_flow_vs_headloss(pipe3_results, 'Pipe 3', 0.5);
writetable(struct2table(pipe3_results), 'pipe3_results.csv');

%---------- Pipe 4 ----------%
fprintf(1, '\nPipe 4\n')
pipe_diameter = 26.04e-3; % 26.04 mm
pipe_length = 1.22; % m
pipe_roughness = 1.5e-6; % 1.5e-3 mm
flow_rate = [8, 16, 24, 32, 40]; % L/min
manometer_reading_mmHg = [0, 0, 0, 3.4, 6.5]; % mmHg
% Run script, get Re, change default values
moody_friction_factor = [0.037, 0.030, 0.026, 0.025, 0.024]; 

pipe4_results = pipe_results(flow_rate, ...
                                       manometer_reading_mmHg, ...
                                       dynamic_viscosity_water, ...
                                       pipe_diameter, ...
                                       pipe_roughness, ...
                                       pipe_length, ...
                                       water_density, ...
                                       moody_friction_factor);

disp(struct2table(pipe4_results));
plot_flow_vs_headloss(pipe4_results, 'Pipe 4', 0.5);
writetable(struct2table(pipe4_results), 'pipe4_results.csv');


%---------- Function Definitions ----------%

function head_loss = head_loss_calc(f, l, u, d)
    % Darcy-Weisbach head loss
    head_loss = (f * l * u^2) / (2 * 9.81 * d);
end

function Re = reynolds_number(density, u, d, dynamic_viscosity)
    Re = (density * u * d) / dynamic_viscosity;
end

function velocity = flow_rate_to_velocity(flow_rate_Lmin, d)
    % Get velocity of the flow
    flow_rate = (flow_rate_Lmin / 60) / 1000; % Convert L/min to m^3/s
    area = (pi / 4) * d^2;
    velocity = flow_rate / area;
end

function friction_factor = haaland_friction_factor(Re, relative_roughness)
    if Re < 4000
        fprintf(1,'Error: Haaland equation is only valid for turbulent flow (Re > 4000)');
        friction_factor = 0;
    else
        friction_factor = (-1.8 * log10( (relative_roughness/3.7)^1.11 + 6.9/Re ))^-2;
    end
end

function friction_factor = colebrook_white_friction_factor(Re, relative_roughness)
% Solve the Colebrook-White equation for the friction factor.

    % Define the function by equating it to zero
    colebrook_eq = @(f) 1./sqrt(f) ...
                              + 2*log10(relative_roughness/3.7 + 2.51./(Re.*sqrt(f)));
    % Use the Haaland equation for an initial guess
    f_guess = haaland_friction_factor(Re, relative_roughness);
    % Finds root on non-linear function
    friction_factor = fzero(colebrook_eq, f_guess); 
end

function err = percent_difference(experimental, theoretical)
    err = abs(experimental - theoretical) / theoretical * 100;
end

function results = pipe_results(flow_rate, ...
                                                manometer_reading_mmHg, ...
                                                dynamic_viscosity, ...
                                                pipe_diameter, ...
                                                pipe_roughness, ...
                                                pipe_length, ...
                                                water_density, ...
                                                f_moody)

    % Generates a table comparing experimental head loss to theoretical
    % Flow rate in pipe as array of values, L/min
    % Manometer reading as array of values, mmHg
    % P1 pressure readings at start of pipe as array of values, mmHg
    % P2 pressure readings at end of pipe as array of values, mmHg
    % Dynamic viscosity, Pa*s
    % Pipe diameter, m
    % Pipe roughness, m
    % Pipe length, m
    % Water density, kg/m^3
    
    fprintf(1, ['Pipe Diameter = %g mm, ' ...
                'Pipe Length = %g mm, ' ...
                'Pipe Roughness = %g mm, ' ...
                'Relative Roughness %.3e\n'], ...
                 pipe_diameter*1000, ...
                 pipe_length*1000, ...
                 pipe_roughness*1000, ...
                 pipe_roughness/pipe_diameter)

    results = []; % initialise
    for i = 1:length(flow_rate)
        velocity = flow_rate_to_velocity(flow_rate(i), pipe_diameter);
        rel_roughness = pipe_roughness / pipe_diameter;
        Re = reynolds_number(water_density, velocity, pipe_diameter, dynamic_viscosity);
        
        f_cole = colebrook_white_friction_factor(Re, rel_roughness);
        f_haaland = haaland_friction_factor(Re, rel_roughness);

        h_l_cole = head_loss_calc(f_cole, pipe_length, velocity, pipe_diameter);
        h_l_haaland = head_loss_calc(f_haaland, pipe_length, velocity, pipe_diameter);
        h_l_moody = head_loss_calc(f_moody(i), pipe_length, velocity, pipe_diameter);
        h_l_exp = manometer_reading_mmHg(i)* 133.322 / (9.81*water_density);

        diff_cole = percent_difference(h_l_exp, h_l_cole);
        diff_haaland = percent_difference(h_l_exp, h_l_haaland);
        diff_moody = percent_difference(h_l_exp, h_l_moody);

        results = [results, struct( ...
            'Flow_L_min', flow_rate(i), ...
            'Mano_mmHg', manometer_reading_mmHg(i), ...
            'Mano_Pa', manometer_reading_mmHg(i) * 133.322, ...
            'Re', Re, ...
            'velocity_m_s', velocity, ...
            'F_Moody', f_moody(i), ...
            'F_Cole', f_cole, ...
            'F_Haaland', f_haaland, ...
            'H_L_Exp_m', h_l_exp, ...
            'H_L_Moody_m', h_l_moody, ...
            'H_L_Cole_m', h_l_cole, ...
            'H_L_Haaland_m', h_l_haaland, ...
            'diff_Moody', diff_moody, ...
            'diff_Cole', diff_cole, ...  
            'diff_Haaland', diff_haaland, ...  
            'H_L_Moody_Pa', h_l_moody * water_density * 9.81, ...
            'H_L_Cole_Pa', h_l_cole * water_density * 9.81, ...
            'H_L_Haaland_Pa', h_l_haaland * water_density * 9.81 ...
            )];
            
    end
end

function flow_figure = plot_flow_vs_headloss(pipe_results, pipe_name, y_max)
    % Creates a plot for flow rate verses the headloss in a pipe.
    % Pipe results is generated with its named function
    % Pipe name is for the title, 'Flow Rate vs. Head Loss for pipe_name'
    % y_max sets the max y value to had consistent scale, but still adjust
    %   for pipes with low head loss.
    flow_figure = figure; % Create new figure window

    % Convert to vectors
    flow = [0, pipe_results.Flow_L_min];
    headloss_exp = [0, pipe_results.H_L_Exp_m];
    headloss_moody = [0, pipe_results.H_L_Moody_m];
    headloss_cole = [0, pipe_results.H_L_Cole_m];
    headloss_haaland = [0, pipe_results.H_L_Haaland_m];
    
    % Plot theoretical as a line and experimental as points
    h_moody = plot(flow, headloss_moody, 'b--', 'DisplayName', 'Moody Chart Head Loss');
    hold on;
    h_cole = plot(flow, headloss_cole, 'g--', 'DisplayName', 'Colebrook-White Head Loss');
    h_haaland = plot(flow, headloss_haaland, 'm--', 'DisplayName', 'Haaland Head Loss');
    h_exp = plot(flow, headloss_exp, 'ro', 'DisplayName', 'Experimental Head Loss');
    hold off;

    title(sprintf('Flow Rate vs. Head Loss for %s', pipe_name));
    
    xlabel('Flow Rate [L/min]'); % X-axis is Flow Rate
    xlim([0 40]);
    xticks([0 8 16 24 32 40]);
    
    ylabel('Headloss [m]'); % Y-axis is Headloss
    ylim([0 y_max]);
    %yticks([0 1 2 3 4]); % Uncomment to set specific ticks
    
    % Grid and minor ticks
    grid on;
    ax = gca;
    ax.YMinorTick = 'on';
    ax.YMinorGrid = 'on';
    ax.XMinorTick = 'on';
    ax.XMinorGrid = 'on';
    
    legend([h_exp, h_moody, h_cole, h_haaland], 'Location', 'northwest');
    
    % Control size and aspect ratio of the plot for consistency
    flow_figure.Units = 'inches';
    flow_figure.Position = [1 1 6 3]; % x-pos, y-pos, x-width, y-width
end