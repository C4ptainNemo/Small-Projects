%% Filter Order Analysis
% Requires the 'Control System Toolbox' and 'Signal Processing Toolbox'.
% This script is designed to analyse the order of filter needed based on
% the cutoff frequency.
% The first section calculates the order needed for a given passband
% frequency, and attenuation in the passband and stopband.
% The second section calculates the max attenuation possible for a chosen
% cutoff frequency.

% ----- Clear Commands ----- %
clc;              % Clear the Command Window
clear;            % Clear variables
close all;        % Close all figures

%% Max Passband Frequency for Nth Order
%----- Filter Specifications -----%
fs = 7800;              % stopband frequency (15625/2), Hz
Amin = -20*log10(1/2^8); % stopband attenuation, approx 48 dB
Amax = 0.1;             % passband attenuation (ripple), dB

% ----- Calculations ----- %
results = []; % results for butt and cheb order per fp
fprintf(1, 'fs = %g Hz, A_min = %.4g dB, A_max = %g dB\n', fs, Amin, Amax)
for j = 100:10:fs % passband frequency, goes from 100->fs in steps of 10 Hz
    fp = j;
    ws = 2*pi*fs;
    wp = 2*pi*fp;

    % Calculate the order and wn for each passband frequency
    [n_butt, wn_butt] = buttord(wp, ws, Amax, Amin, 's'); 
    [n_cheb, wn_cheb] = cheb1ord(wp, ws, Amax, Amin, 's');

    results = [results, struct( ...
        'fp', fp, ...
        'n_butt', n_butt, ...
        'n_cheb', n_cheb ...
        )];
end
%disp(struct2table(results))

% ----- Find max fp for a given filter order ----- %
butt_ord = [];
cheb_ord = [];
for j = 2:length(results) % start at 2 since the previous index is checked
% Check if the current and previous orders are the same
% If not, append the previous order and its freq to the struct
    if results(j).n_butt ~= results(j-1).n_butt
        butt_ord = [butt_ord, struct( ...
        'n_butt', results(j-1).n_butt, ...
        'fp_max', results(j-1).fp)];
    end

    if results(j).n_cheb ~= results(j-1).n_cheb
        cheb_ord = [cheb_ord, struct( ...
        'n_cheb', results(j-1).n_cheb, ...
        'fp_max', results(j-1).fp)];
    end
end
%disp(struct2table(butt_min_ord)) % Displays the data as a table
%disp(struct2table(cheb_min_ord)) % Displays the data as a table

% ---- Plot Filter Orders for Max fp ----- %
x1 = [butt_ord.fp_max]; % Turn struct into vectors for plotting
y1 = [butt_ord.n_butt];
x2 = [cheb_ord.fp_max];
y2 = [cheb_ord.n_cheb];

fig1 = figure; % Create new figure to plot on
scatter(x1, y1, 'magenta', '*');    % Butterworth Plot
hold on;
scatter(x2, y2, 'black', 'o');      % Chebyshev Plot

xlabel('Frequency [Hz]');
ylabel('Filter Order');
legend('Butterworth', 'Chebyshev', 'Location', 'northwest');
title(sprintf('Maximum Cutoff Frequency Per Filter Order with %g dB Ripple', Amax));
grid on;
set(gca, 'YScale', 'log'); % Sets the y-axis to log scale, useful for
%                             showing the full scale, since the filter 
%                             order grows exponentially as fp approached fs
ax = gca;               % get current axis
ax.YMinorTick = 'on';   % Turn on minor ticks
ax.YMinorGrid = 'on';   % Turn on minor grid lines
%ax.YLim = [0, 20];      % Set y limit, max ends up being in the 200's
ax.XMinorTick = 'on';   % Turn on minor ticks
ax.XMinorGrid = 'on';   % Turn on minor grid lines

% Control size and aspect ratio of the plot for consistency
% Make sure the x-y ratios are the same for both functions
fig1.Units = 'inches';
fig1.Position = [1 1 12 4];  % Plot size, x-pos, y-pos, x-width, y-width

hold off;

% ----- Plot Transfer Functions ----- %
% For each fp in in butt/cheb_min_ord loop through and plot the 
%   transfer function

fig2 = figure; % Create new figure to plot on
max_ord = 20;
for j = 1:max_ord
    fp = cheb_ord(j).fp_max;      % passband frequency, Hz
    ws = 2*pi*fs;
    wp = 2*pi*fp;
    
    [~, wn_cheb] = cheb1ord(wp, ws, Amax, Amin, 's');
    [b, a] = cheby1(cheb_ord(j).n_cheb, Amax, wn_cheb,'low','s');
    H = tf(b,a);
    h = bodeplot(H);
    hold on;
end

setoptions(h, ...
    'YLim', {[-60, 5], [-1800, 10]}, ...
    'XLim', [100, 10000], ...
    'FreqUnits','Hz', ...
    'PhaseVisible','on');
title(sprintf(['Frequency Response of Order 1 to %g Chebyshev Filters at Maximum ' ...
    'Passband Frequency with %g dB Ripple'], max_ord, Amax));
grid on;
ax = gca;               % get current axis
ax.XMinorTick = 'on';   % Turn on minor ticks
ax.XMinorGrid = 'on';   % Turn on minor grid lines
ax.YMinorTick = 'on';   % Turn on minor ticks
ax.YMinorGrid = 'on';   % Turn on minor grid lines

% Control size and aspect ratio of the plot for consistency
% Make sure the x-y ratios are the same for both functions
fig2.Units = 'inches';
fig2.Position = [1 1 12 8];  % Plot size, x-pos, y-pos, x-width, y-width

hold off;

%% Max Attenuation For Given Order and Passband Frequency
% Set the desired passband frequency fp.
% This will then go through each order of filter and find the max
%   attenuation possible before reaching the stopband frequency.
% Do for both Butterworth and Cheyshev. Was easiest to just copy paste code.

fp = 3400; % Hz, Set this based off previous sections analysis
wp = 2*pi*fp; % rad/s
max_order = 20; % Maximum order of filter to check

butt_attenuation_results(max_order) = struct('Order', [], 'Amin', []);
% Butterworth
for i = 1:max_order
    order = i;
    for j = 1:1000 % Will check up to 1000 dB of attenuation
    % Iterate through attenuation, calculate the order needed.
    % Once n_butt exceeds order exit loop.
        Amin = j;
        [n_butt, ~] = buttord(wp, ws, Amax, Amin, 's');
        if n_butt > order % Once calculated order exceeds order
            butt_attenuation_results(i).Order = order;
            butt_attenuation_results(i).Amin = Amin - 1;
            break
        end
    end
end
%disp(struct2table(butt_attenuation_results))

% Chebyshev
cheb_attenuation_results(max_order) = struct('Order', [], 'Amin', []);
for i = 1:max_order
    order = i;
    for j = 1:1000 % Will check up to 200 dB of attenuation
    % Iterate through attenuation, calculate the order needed.
    % Once n_butt exceeds order exit loop.
        Amin = j;
        [n_cheb, ~] = cheb1ord(wp, ws, Amax, Amin, 's');
        if n_cheb > order % Once calculated order exceeds order
            cheb_attenuation_results(i).Order = order;
            cheb_attenuation_results(i).Amin = Amin - 1;
            break
        end
    end
end
%disp(struct2table(cheb_attenuation_results))

fig3 = figure; % Create new figure to plot on
x1 = [butt_attenuation_results.Order]; % Turn struct into vectors for plotting
y1 = [butt_attenuation_results.Amin];
x2 = [cheb_attenuation_results.Order];
y2 = [cheb_attenuation_results.Amin];
scatter(x1, y1, 'magenta', '*');   % Butterworth Plot
hold on;
scatter(x2, y2, 'black', 'o');   % Chebyshev Plot

xlabel('Filter Order');
ylabel('Maximum Attenuation [dB]');
legend('Butterworth', 'Chebyshev', 'Location', 'northwest');
title(sprintf('Max Attenuation for Nth Order Filter with %g dB Ripple', Amax));
grid on;
ax = gca;               % get current axis
ax.YMinorTick = 'on';   % Turn on minor ticks
ax.YMinorGrid = 'on';   % Turn on minor grid lines
ax.XMinorTick = 'on';   % Turn on minor ticks
ax.XMinorGrid = 'on';   % Turn on minor grid lines

% Control size and aspect ratio of the plot for consistency
% Make sure the x-y ratios are the same for both functions
fig3.Units = 'inches';
fig3.Position = [1 1 12 4];  % Plot size, x-pos, y-pos, x-width, y-width

hold off;