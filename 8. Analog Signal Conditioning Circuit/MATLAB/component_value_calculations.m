%% Component Value Calculations
% Requires the 'Control System Toolbox' and 'Signal Processing Toolbox'.
% This script is used to calculate the component values for an analog lowpass 
%   filter with a Chebyshev response and VCVS topology.
% Should first run through Filter_order_analysis to get help work out filter specs.
% I recommend hitting publish to more easily see the large tables generated.

% ----- Table Notes ----- %
% The 'component_calculator_sallen_key' functions will output a table of
% tuned component values.
% k: Tuning iteration.
% wn_ack: Natural frequency of stage from rounded component values.
% Q_ack: Q-factor of stage from rounded component values.
% wn_Error: Percentage error between input wn and wn_ack.
% Q_Error: Percentage error between input Q and Q_ack.
% R1, R2, C1, C2: Ideal component values for the stage.
% R1_Rounded -> C2_rounded: Closest E12 components to the ideal ones.
% R1_Error -> C2_Error: Percentage error between ideal and rounded components.
% Avg_Error: Average of R1_Error, R2_Error and C2_Error. C1's error is zero.

% ----- Clear Commands ----- %
clc;              % Clear the Command Window
clear;            % Clear variables
close all;        % Close all figures

% ----- Filter Specifications ----- %
fs = 7800;      % stopband frequency (15625/2), Hz
fp = 3400;      % passband frequency, Hz
ws = 2*pi*fs;   % stopband frequency, rad/s
wp = 2*pi*fp;   % passband frequency, rad/s

C_base = 2.2e-9; % Set value for C1 of each stage, rule of thumb is 10/fp uF
                 % 10/3400 uF ~ 2.94 nF

Amin = -20*log10(1/2^8);  % ~ 48 dB
Amax = 0.1; % passband attenuation (ripple), dB

fprintf(1,['\nThree Stage Chebyshev Response, ' ...
            'Sallen-Key Topology Lowpass Filter\n' ...
            'Passband Frequency: %g Hz\n' ...
            'Stopband Frequency: %g Hz\n' ...
            'Maximum Passband Gain: %g dB\n' ...
            'Minimum Stopband Attenuation: %.4g dB\n\n'],...
    fp,fs,Amax,Amin);

% ----- Calculate Chebyshev Order & Transfer Function ----- %
[n_cheb, wn_cheb] = cheb1ord (wp, ws, Amax, Amin, 's'); % Order and natural freq
[b, a] = cheby1(n_cheb, Amax, wn_cheb,'low','s'); % Transfer func coefficients
H = tf(b,a); % Transfer Function
[z, p, k] = tf2zpk (b, a); % Get the poles
[magnitude_at_ws, ~, ~] = bode(H, ws); % Takes H and ws and returns magnitude
fprintf(1, 'Chebyshev Filter Order: %g\n', n_cheb); % Sanity check it's 6th Order
fprintf(1, 'Transfer Function Attenutation at %g Hz = %.4g dB\n\n', ...
    fs, -20*log10(magnitude_at_ws)); % Convert magnitude to dB

%%%%% stage A %%%%%
% Use poles to get 2nd order transfer function coefficients.
% Use the coefficients to get the Q-factor and natural frequency.
% Input those along with a base capacitor values to get a table of tuned components.
a1A = -p(1)-p(2);
a2A = p(1)*p(2);
wnA = sqrt(a2A);
QA = wnA/a1A;
fprintf(1,'Stage A: wn = %g rad/s, fn = %g Hz, Q = %g\n', ...
           wnA, wnA/(2*pi), QA);
stageAresults = component_calculator_vcvs(QA, wnA, C_base, true, false);

%%%%% stage B %%%%
a1B = -p(3)-p(4);
a2B = p(3)*p(4);
wnB = sqrt(a2B);
QB = wnB/a1B;
fprintf(1,'Stage B: wn = %g rad/s, fn = %g Hz, Q = %g\n', ...
           wnB, wnB/(2*pi), QB);
stageBresults = component_calculator_vcvs(QB, wnB, C_base, true, false);

%%%%% stage C %%%%
a1C = -p(5)-p(6);
a2C = p(5)*p(6);
wnC = sqrt(a2C);
QC = wnC/a1C;
fprintf(1,'Stage C: wn = %g rad/s, fn = %g Hz, Q = %g\n', ...
           wnC, wnC/(2*pi), QC);
stageCresults = component_calculator_vcvs(QC, wnC, C_base, true, false);

%% ----- Compare Transfer Functions -----
% Takes the real rounded component values of each stage and calculates the
% transfer function, then multiplies all three to get the final function. 
% The theoretical one calculated at the start is plotted against the real 
% one.

% Select the desired line for each stage (input k+1)
stageAline = 13; % 13 for actual value, 12 for high error plot
stageBline = 24; % 24 for actual value, 16 for high error plot
stageCline = 17; % 17 for actual value, 13 for high error plot

% This outputs a table showing the three selected lines
a = stageAresults(stageAline); % Extract desired line
b = stageBresults(stageBline);
c = stageCresults(stageCline);
fields = fieldnames(a); % Get the fields for each column
combined_stage_lines = struct();
for i = 1:length(fields)
    field_name = fields{i}; % Get the actual field name string
    combined_stage_lines(1).(field_name) = a.(field_name);
    combined_stage_lines(2).(field_name) = b.(field_name);
    combined_stage_lines(3).(field_name) = c.(field_name);
end
fprintf(1, 'Table showing the stage lines selected\n')
disp(struct2table(combined_stage_lines));


% Calculate the transfer function for each stage
H_stageA = transfer_function_vcvs(1, ...
    stageAresults(stageAline).R1_rnd, stageAresults(stageAline).R2_rnd,...
    stageAresults(stageAline).C1_rnd, stageAresults(stageAline).C2_rnd);

H_stageB = transfer_function_vcvs(1, ...
    stageBresults(stageBline).R1_rnd, stageBresults(stageBline).R2_rnd,...
    stageBresults(stageBline).C1_rnd, stageBresults(stageBline).C2_rnd);

H_stageC = transfer_function_vcvs(1, ...
    stageCresults(stageCline).R1_rnd, stageCresults(stageCline).R2_rnd,...
    stageCresults(stageCline).C1_rnd, stageCresults(stageCline).C2_rnd);

% Combine the transfer functions
H_Calculated = H_stageA * H_stageB * H_stageC;

% Gives the magnitude at the natural frequencies
[magA, ~, ~] = bode(H_stageA, wnA);
magA_db = 20*log10(magA);
[magB, ~, ~] = bode(H_stageB, wnB);
magB_db = 20*log10(magB);
[magC, ~, ~] = bode(H_stageC, wnC);
magC_db = 20*log10(magC);
fprintf(1, ['Gain at natural frequencies: Stage A = %g dB, Stage B = %g dB, ' ...
        'Stage C = %g dB\n'], ...
        magA_db, magB_db, magC_db);

% Plot theoretical vs calculated transfer functions
fig1 = figure;
fig1.Units = 'inches';
fig1.Position = [1 1 12 6];  % x, y, width, height (in inches)
h = bodeplot(H, 'g-', H_Calculated, 'r--');
title('Theoretical Vs Calculated Transfer Functions');
setoptions(h, ...
    'FreqUnits','Hz', ...
    'PhaseVisible','on', ...
    'XLim', [10, 10000], ...
    'YLim', {[-60, 10], [-540, 10]});
grid on;
legend('Theoretical', 'Calculated', 'Location', 'southwest');


%% ----- Functions -----
% Contains small functions used above
% component_calculator_vcvs is large so gets own file.

function transfer_function = transfer_function_vcvs(K, R1, R2, C1, C2)
% Takes component values for a vcvs key lowpass filter and returns 
% the transfer function.
    num = [0 0 K/(R1*R2*C1*C2)];    % numerator
    den = [1 ...                    % denominator
           1/(R1*C2) + 1/(R2*C2) + (1-K)/(R2*C1) ...
           1/(R1*R2*C1*C2)]; 
    transfer_function = tf(num,den);
end