%% EGB211 Computer Lab Task
%
% Written by Adam Telfer, n12068373.
%
% All units are in SI.
%
% 'variable'_0 indicates initial state of variable.
% 'variable'_dt indicates the first time derivative.
% 'variable'_dt2 indicates the second time derivative.
%
% Plotting for each question is split into it's own block for visual
% clarity.

clear; close all; clc;

%% Finite Difference Method Function
% Used for questions 3, 4 and 7.

function [x, y] = finite_diff_q3_4_7(g, m, k, dt, N_steps, ...
                                     L_0, r_0, r_dt_0, theta_0, theta_dt_0)
    % Initialise arrays
    x = zeros(N_steps, 1);
    x_dt = zeros(N_steps, 1);
    x_dt2 = zeros(N_steps, 1);
    
    y = zeros(N_steps, 1);
    y_dt = zeros(N_steps, 1);
    y_dt2 = zeros(N_steps, 1);
    
    % Initialise x and y positions
    x(1) = r_0 * sin(theta_0);
    x(2) = x(1) + dt * ...
        (r_dt_0 * sin(theta_0) + r_0 * theta_dt_0 * cos(theta_0));
    
    y(1) = r_0 * cos(theta_0);
    y(2) = y(1) + dt * ...
        (r_dt_0 * cos(theta_0) - r_0 * theta_dt_0 * sin(theta_0));
    
    % Run finite difference simulation.
    for j = 2:N_steps - 1
        % Assign values for previous and current positions.
        x_PRE = x(j-1);
        x_CUR = x(j);
        y_PRE = y(j-1);
        y_CUR = y(j);
        
        % Calculate current radial length and spring force.
        r = sqrt(x_CUR^2 + y_CUR^2);
        F_spring = k * (r - L_0);
    
        % Calculate spring force components.
        % sin(theta) = x/r and cos(theta) = y/r.
        % The sign of the force comes from the sign of x_CUR or y_CUR.
        F_x = F_spring * (x_CUR/r);
        F_y = F_spring * (y_CUR/r);
        
        % Functions for the equations of motion in the x and y directions,
        % equated to zero.
        f_x = @(x_NEX) ...
            (m*((x_NEX - 2*x_CUR + x_PRE)/(dt^2)) + F_x);
        f_y = @(y_NEX) ...
            (m*((y_NEX - 2*y_CUR + y_PRE)/(dt^2)) + F_y - m*g);
    
        % Find solutions to the equations
        x_NEX = fzero(f_x, x_CUR);
        y_NEX = fzero(f_y, y_CUR);
    
        % Store results for position, velocity, acceleration.
        x(j+1) = x_NEX;
        x_dt(j) = (x_NEX - x_PRE) / (2*dt);
        x_dt2(j) = (x_NEX - 2 * x_CUR + x_PRE) / (dt^2);
        
        y(j+1) = y_NEX;
        y_dt(j) = (y_NEX - y_PRE) / (2*dt);
        y_dt2(j) = (y_NEX - 2 * y_CUR + y_PRE) / (dt^2);
    end
end

%% Q3 - Finite Difference

% Constants
g = 9.8; % m/s^2, Gravitational acceleration.
m = 0.1; % kg, Mass of weight.
k = 20; % N/m, Elastic coefficient of spring.

% Initial State
L_0 = 20; % m, Unstretched length of spring.

r_0 = 22; % m, Initial length of spring.
r_dt_0 = 0.1; % m/s, Initial radial velocity.

theta_0 = 0.05; % rad, Initial angle.
theta_dt_0 = 0.01; % rad/s, Initial angular velocity.

% Simulation Parameters
dt = 0.1; % s, Time step.
t_end = 10; % s, Total time of simulation.
N_steps = round(t_end / dt, 0); % Number of time steps
    
% Run simulation
[x_Q3, y_Q3] = finite_diff_q3_4_7(g, m, k, dt, N_steps, ...
                                  L_0, r_0, r_dt_0, theta_0, theta_dt_0);
    
%% Plotting for Q3

fig_q3 = figure;

% Plot path, uses the x and y vectors stored from the first iteration.
plot(x_Q3, y_Q3, 'DisplayName', 'Path');
hold on;
% Marker for initial position
plot(x_Q3(1),y_Q3(1), 'ro', 'MarkerSize', 10, ...
    'DisplayName', 'Initial Position');
title('Path of Pendulum');
xlabel('x [m]');
ylabel('y [m]');
grid on;
set(gca, 'YDir', 'reverse'); % Flips the y-axis
legend('Location', 'Southeast');
axis equal;

% Control the size of the figure
set(fig_q3, 'Units', 'inches', 'Position', [1, 1, 16, 8]);

% Save the figure as a .tiff file
exportgraphics(gcf, 'plot_Q3.tif', 'Resolution', 600);

%% Q4 - Convergence Study

% Runs a convergence study on the simulation.
% Outputs table of results in the command window and also saves a .xlsx
% file to make it easy to copy and paste.

% Constants
g = 9.8; % m/s^2, Gravitational acceleration.
m = 0.1; % kg, Mass of weight.
k = 20; % N/m, Elastic coefficient of spring.

% Initial State
L_0 = 20; % m, Unstretched length of spring.

r_0 = 22; % m, Initial length of spring.
r_dt_0 = 0.1; % m/s, Initial radial velocity.

theta_0 = 0.05; % rad, Initial angle.
theta_dt_0 = 0.01; % rad/s, Initial angular velocity.

% Simulation Parameters
dt = 0.1; % s, Time step.
t_end = 10; % s, Total time of simulation.

convergence_cutoff = 0.005; % m, cutoff to stop convergence study when the 
% change in final position between the previous and current iteration are 
% below this value.
convergence_max_iterations = 100; % Stops infinite loop if convergence 
% takes too long or doesn't occur if a while loop was used instead of
% a for loop.

% Outer loop (i variable) is for convergence study.
% Inner loop (j variable) is for the finite difference simulation.
for i = 1:convergence_max_iterations
    tic % Record time each iteration takes.
    N_steps = round(t_end / dt, 0); % Number of time steps
    
    % Run simulation
    [x, y] = finite_diff_q3_4_7(g, m, k, dt, N_steps, ...
                                L_0, r_0, r_dt_0, theta_0, theta_dt_0);
    
    simulation_elapsed_time = toc; % Store time iteration took

    % Store the time step, number of time steps, time iteration took,
    % last x and y position, change in final position between current and
    % previous iteration.
    convergence_results.dt(i) = dt;
    convergence_results.N_steps(i) = N_steps;
    convergence_results.run_time(i) = simulation_elapsed_time;
    convergence_results.x(i) = x(N_steps);
    convergence_results.y(i) = y(N_steps);
    if i == 1
        % There will be no change in the first iteration.
        convergence_results.change(i) = NaN;
    else
        % Calculate the change in the x and y positions.
        dx = convergence_results.x(i-1) - convergence_results.x(i);
        dy = convergence_results.y(i-1) - convergence_results.y(i);
        ds = sqrt(dx^2 + dy^2);
        convergence_results.change(i) = ds;
    end

    % If the change in final position between current and previous 
    % iteration is less than the cutoff value, exit the loop.
    if convergence_results.change(i) < convergence_cutoff
        break
    end

    % Half the time step for the next iteration.
    dt = dt / 2;
end

x_Q4 = x;
y_Q4 = y;

% Create a table to display the final results.
x_title = sprintf('x(t=%.2f) [m]', t_end);
y_title = sprintf('y(t=%.2f) [m]', t_end);
convergence_results_table = table( ...
    convergence_results.dt', ...
    convergence_results.N_steps', ...
    convergence_results.run_time', ...
    convergence_results.x', ...
    convergence_results.y', ...
    convergence_results.change', ...
    'VariableNames', {'dt [s]', 'N Steps', 'Run Time [s]', x_title, ...
                       y_title, 'Change [m]'});

% Display the value of parameters and the table
fprintf('Simulation Length: %.2f [s]\n', t_end);
fprintf('Convergence Cutoff: %.3f [m]\n', convergence_cutoff);

fprintf('Initial Radius: %.2f [m]\n', r_0);
fprintf('Initial Radial Velocity: %.2f [m/s]\n', r_dt_0);

fprintf('Initial Angle: %.2f [rad]\n', theta_0);
fprintf('Initial Angular Velocity: %.2f [rad/s]\n', theta_dt_0);

fprintf('Spring Constant: %.2f [N/m]\n', k);
fprintf('Natural Length of Spring: %.2f [m]\n', L_0);
fprintf('Mass: %.2f [kg]\n', m);
fprintf('Gravitational Acceleration: %.2f [m/s^2]\n', g);

fprintf('\n')
disp(convergence_results_table);

% Save the convergence results table to an Excel file.
writetable(convergence_results_table, 'convergence_table.xlsx');

%% Plotting for Q4

% Create time vector for plotting.
time_vector = (0:N_steps-1) * dt;

fig_q4 = figure;
sgtitle('Position of Elastic Pendulum')

% Uses the x and y vectors of the last iteration
% Subplot for x position
subplot(1, 3, 1);
plot(time_vector, x_Q4);
title('x-Position Over Time');
xlabel('Time [s]');
ylabel('x [m]');
grid on;

% Subplot for y position
subplot(1, 3, 2);
plot(time_vector, y_Q4);
title('y-Position Over Time');
xlabel('Time [s]');
ylabel('y [m]');
grid on;

% Subplot for path
subplot(1, 3, 3);
% Plot path
plot(x_Q4, y_Q4, 'DisplayName', 'Path of Mass');
hold on;
% Marker for initial position
plot(x(1),y(1), 'ro', 'MarkerSize', 10, 'DisplayName', 'Initial Position');
title('Path');
xlabel('x [m]');
ylabel('y [m]');
grid on;
set(gca, 'YDir', 'reverse'); % Flips the y-axis
legend('Location', 'Southeast');
axis equal;

% Control the size of the figure
set(fig_q4, 'Units', 'inches', 'Position', [1, 1, 16, 8]);

% Save the figure as a .tiff file
exportgraphics(gcf, 'plot_Q4.tif', 'Resolution', 600);

%% Q7 - Model Validation

% The analytical model and plotting have been made into functions so that
% it is easier to run different conditions and parameters.

% Analytical Model of Elastic Pendulum
function [x, y] = analytical_model(g, m, k, dt, t_end, ...
                                   L_0, r_0, r_dt_0, theta_0, theta_dt_0)
    % Analytical model of the elastic pendulum, done using polar 
    % coordinates. Can be compared to the numerical solution to validate 
    % it.
    %
    % Only valid for small angles, since small angle approximation was used
    % to turn a non-linear differential equation ionto a linear one to be 
    % able to solve it.

    % Create a time vector for the analytical solution.
    time_vector = 0:dt:t_end;
    
    % Derived analytical equation for theta direction.
    w_n_theta = sqrt(g / L_0);
    A_theta = sqrt(theta_0^2 + (theta_dt_0 / w_n_theta)^2);
    phi_theta = atan(theta_dt_0 / (theta_0 * w_n_theta));
    theta_analytical = A_theta * cos(w_n_theta * time_vector - phi_theta);
    
    % Derived analytical equation for the radial direction.
    w_n_radial = sqrt(k / m);
    C_radial = L_0 + m * g / k;
    A_radial = sqrt((r_0 - C_radial)^2 + (r_dt_0 / w_n_radial)^2);
    phi_radial = atan(r_dt_0 / ((r_0 - C_radial) * w_n_radial));
    radial_analytical = A_radial * ...
        cos(w_n_radial * time_vector - phi_radial) + C_radial;
    
    % To plot the analytical solution convert the polar coordinates into
    % rectangular.
    x = radial_analytical .* sin(theta_analytical);
    y = radial_analytical .* cos(theta_analytical);
end

% Constants
g = 9.8; % m/s^2, Gravitational acceleration.

% Initial State
L_0 = 20; % m, Unstretched length of spring.

r_0 = 22; % m, Initial length of spring.
r_dt_0 = 0.1; % m/s, Initial radial velocity.

theta_0 = 0.05; % rad, Initial angle.
theta_dt_0 = 0.01; % rad/s, Initial angular velocity.

% Simulation Parameters
dt = 0.003; % s, Time step.
  
% Simulation 1
m = 0.1; % kg, Mass of weight.
k = 20; % N/m, Elastic coefficient of spring.
t_end = 10; % s, Total time of simulation.
N_steps = round(t_end / dt, 0); % Number of time steps

[x_Q7_1_numerical, y_Q7_1_numerical] = ...
    finite_diff_q3_4_7(g, m, k, dt, N_steps, ...
                       L_0, r_0, r_dt_0, theta_0, theta_dt_0);

[x_Q7_1_analytical, y_Q7_1_analytical] = ...
    analytical_model(g, m, k, dt, t_end, ...
                     L_0, r_0, r_dt_0, theta_0, theta_dt_0);

plot_Q7(x_Q7_1_numerical, y_Q7_1_numerical, ...
        x_Q7_1_analytical, y_Q7_1_analytical, dt, N_steps, t_end, 1);

% Simulation 2
m = 10; % kg, Mass of weight.
k = 10; % N/m, Elastic coefficient of spring.
t_end = 60; % s, Total time of simulation.
N_steps = round(t_end / dt, 0); % Number of time steps

[x_Q7_2_numerical, y_Q7_2_numerical] = ...
    finite_diff_q3_4_7(g, m, k, dt, N_steps, ...
                       L_0, r_0, r_dt_0, theta_0, theta_dt_0);

[x_Q7_2_analytical, y_Q7_2_analytical] = ...
    analytical_model(g, m, k, dt, t_end, ...
                     L_0, r_0, r_dt_0, theta_0, theta_dt_0);

plot_Q7(x_Q7_2_numerical, y_Q7_2_numerical, ...
        x_Q7_2_analytical, y_Q7_2_analytical, dt, N_steps, t_end, 2);

%% Plotting for Q7
function fig_q7 = plot_Q7(x_numerical, y_numerical, ...
    x_analytical, y_analytical, dt, N_steps, t_end, simulation_number)
    filename = sprintf('plot_Q7_%i.tif', simulation_number);
    
    % Create time vector for plotting.
    time_vector_numerical = (0:N_steps-1) * dt;
    time_vector_analytical = 0:dt:t_end;
    
    fig_q7 = figure;
    sgtitle('Position of Elastic Pendulum')
    
    % Subplot for x position
    subplot(1, 3, 1);
    % Numerical Plot
    plot(time_vector_numerical, x_numerical);
    hold on;
    % Analytical Plot
    plot(time_vector_analytical, x_analytical, 'r--')
    hold off;
    title('x-Position Over Time');
    xlabel('Time [s]');
    ylabel('x [m]');
    grid on;
    
    % Subplot for y position
    subplot(1, 3, 2);
    % Numerical Plot
    plot(time_vector_numerical, y_numerical);
    hold on;
    % Analytical Plot
    plot(time_vector_analytical, y_analytical, 'r--');
    hold off;
    title('y-Position Over Time');
    xlabel('Time [s]');
    ylabel('y [m]');
    grid on;
    
    % Subplot for path
    subplot(1, 3, 3);
    % Numerical Plot
    plot(x_numerical, y_numerical, 'DisplayName', 'Numerical');
    hold on;
    % Analytical Plot
    plot(x_analytical, y_analytical, 'r--', 'DisplayName', 'Analytical')
    % Marker for initial position
    plot(x_numerical(1),y_numerical(1), 'go', 'MarkerSize', 10, ...
        'DisplayName', 'Initial Position');
    hold off;
    title('Path');
    xlabel('x [m]');
    ylabel('y [m]');
    grid on;
    set(gca, 'YDir', 'reverse'); % Flips the y-axis
    legend('Location', 'Southeast');
    axis equal;
    
    % Control the size of the figure
    set(fig_q7, 'Units', 'inches', 'Position', [1, 1, 16, 8]);
    
    % Save the figure as a .tif file
    exportgraphics(gcf, filename, 'Resolution', 600);
end

%% Q8 - Bungee Jumper

% Simulates bungee jumper. Can change parameters to see what happens, like
% setting m = 80 will cause them to hit the ground.

% Constants
g = 9.8; % m/s^2, gravitational acceleration
m = 60; % kg, mass of bungee jumper
k = 20; % N*m, elastic coefficient of bungee cord

L_0 = 20; % m, Unstretched length of bungee cord

% Initial Conditions
x_0 = 2; % m, Initial position of bungee jumper relative to fixed end
y_0 = 0; % m, Initial position of bungee jumper relative to fixed end

h_0 = 100; % m, Height of fixed end above ground

% Simulation parameters
dt = 0.003; % s, time step for bungee jump simulation
t_end = 120; % s, total time of simulation
N_steps = round(t_end / dt, 0); % Number of time steps

% Initialise arrays
x = zeros(N_steps, 1);
x_dt = zeros(N_steps, 1);
x_dt2 = zeros(N_steps, 1);

y = zeros(N_steps, 1);
y_dt = zeros(N_steps, 1);
y_dt2 = zeros(N_steps, 1);

velocity = zeros(N_steps, 1);
acceleration = zeros(N_steps, 1);
T_cord_array = zeros(N_steps, 1);

% Initialise x and y positions
% This is based on the assumption that the bungee isn't under tension yet.
% If r > L_0 then it won't make much difference anyway since tension will 
% turn on the next time step.
x(1) = x_0;
x(2) = x(1) + 0;

y(1) = y_0;
y(2) = y(1) + 0.5 * g * dt^2;

% Run finite difference simulation.
tic % Record time each iteration takes.
for i = 2:N_steps - 1
    % Assign values for previous and current positions.
    x_PRE = x(i-1);
    x_CUR = x(i);
    y_PRE = y(i-1);
    y_CUR = y(i);
    
    % Calculate current radial length and sigmoid equation.
    r = sqrt(x_CUR^2 + y_CUR^2);
    sigmoid = 1 / (1 + exp(-(r - L_0) ) );
    
    % The velocity in the y-direction is needed to calculate the bungee
    % cord tension. So the equation of motion in the y-direction is used to
    % find the change in the y-position, which then gives the velocity and
    % thus the equation of motion in the x-direction can be solved.

    % Equation of motion and the solution in the y-direction.
    f_y = @(y_NEX) (m*((y_NEX - 2*y_CUR + y_PRE)/(dt^2)) +...
     (k * (r - L_0) + 5*( (y_NEX - y_PRE)/(2*dt)) )...
     * sigmoid * (y_CUR / r) - m * g);
    y_NEX = fzero(f_y, y_CUR);

    % Store position, velocity and acceleration in the y-direction.
    y(i+1) = y_NEX;
    y_dt(i) = (y_NEX - y_PRE) / (2*dt);
    y_dt2(i) = (y_NEX - 2 * y_CUR + y_PRE) / (dt^2);

    % Equation of motion and the solution in the x-direction.
    f_x = @(x_NEX) (m*((x_NEX - 2*x_CUR + x_PRE)/(dt^2)) +...
     (k * (r - L_0) + 5*y_dt(i)) * sigmoid * (x_CUR / r));
    x_NEX = fzero(f_x, x_CUR);
    
    % Store position, velocity and acceleration in the x-direction.
    x(i+1) = x_NEX;
    x_dt(i) = (x_NEX - x_PRE) / (2*dt);
    x_dt2(i) = (x_NEX - 2 * x_CUR + x_PRE) / (dt^2);
    
    % Calculate and store magnitude of velocity and acceleration.
    velocity(i) = sqrt(x_dt(i)^2 + y_dt(i)^2);
    acceleration(i) = sqrt(x_dt2(i)^2 + y_dt2(i)^2);

    % Calculate and store the tension in the bungee cord.
    T_cord = (k * (r - L_0) + 5*y_dt(i)) * sigmoid;
    T_cord_array(i) = T_cord;
end

%% Plotting for Q8

% Create time vector for plotting.
time_vector = (0:N_steps-1) * dt;

% Height, Velocity, and Acceleration Plots
fig_q8 = figure;
sgtitle('Height, Velocity, Acceleration, and Path of Bungee Jumper');

% Subplot for height above ground
subplot(2, 2, 1);
% The -y + h_0 is due to the y values being relative to the fixed point, 
% so this changes it relative to the ground.
plot(time_vector, -y + h_0);
title('Height Above Ground');
xlabel('Time [s]');
ylabel('Height [m]');
grid on;
ylim([0 max(-y + h_0)]); % Set y-axis limits

% Subplot for velocity
subplot(2, 2, 2);
plot(time_vector(2:end-1), velocity(2:end-1));
title('Velocity');
xlabel('Time [s]');
ylabel('Velocity [m/s]');
grid on;
% Set y-axis limits, rounds up to the nearest 5
ylim([0 ceil(max(velocity(2:end-1))/5)*5]); 

% Subplot for acceleration
subplot(2, 2, 3);
plot(time_vector(2:end-1), acceleration(2:end-1));
title('Acceleration');
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]'); 
grid on;
% Set y-axis limits, rounds up to the nearest 5
ylim([0 ceil(max(acceleration(2:end-1))/5)*5]);

% Subplot for path in x and y
subplot(2, 2, 4);
plot(x, y);
title('Path of the Bungee Jumper');
xlabel('x Position [m]');
ylabel('y Position [m]');
grid on;
set(gca, 'YDir', 'reverse');
xlim([min(x)-0.5 max(x)+0.5]); % Set x-axis limits
ylim([min(y) h_0]); % Set y-axis limits

% Control the size of the figure
set(fig_q8, 'Units', 'inches', 'Position', [1, 1, 16, 8]);

% Save the figure as a .tiff file
exportgraphics(gcf, 'plot_Q8.tif', 'Resolution', 600);