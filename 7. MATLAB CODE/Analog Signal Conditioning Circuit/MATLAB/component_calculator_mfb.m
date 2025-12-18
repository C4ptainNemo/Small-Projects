function results = component_calculator_mfb(Q, wn, C_base, print_table, print_ltspice)
% Function for component tuning for a Multiple Feed Back filter.
% Q is the Q-factor for the filter
% wn is the natural frequency of the filter, rad/s
% C_base is the base capacitor value, fro which the other component values
%   are derived from, Farads
% print_table is a Boolean that will display the results from the function,
%   default should be true
% print_ltspice is a Boolean that will display a spice directive input to
%   simulate all the sets of components. Check the section below for more
%   info.

    results = []; % Initialize output
    for k = 0:24
        m = 10^(k / 24);
        C2 = m*C_base;
        R = 1/ (3*m*C_base*wn*Q);
        C1 = C_base*m*(9*Q^2);

        % Round to E12
        R_rounded = round_to_E12(R);
        C1_rounded = round_to_E12(C1);
        C2_rounded = round_to_E12(C2);

        wn_ack = 1 / sqrt(R_rounded*R_rounded*C1_rounded*C2_rounded);
        Q_ack = wn/((3/R_rounded)*(1/C1_rounded));
        wn_percentage_error = percent_error(wn, wn_ack);
        Q_percentage_error = percent_error(Q, Q_ack);

        R_pecentage_error = percent_error(R, R_rounded);
        C1_pecentage_error = percent_error(C1, C1_rounded);
        C2_pecentage_error = percent_error(C2, C2_rounded);
        Average_error = (R_pecentage_error + C1_pecentage_error + C2_pecentage_error) / 3;

        % Store results
        results = [results; struct(...
            'k', k, ...
            'wn_ack', wn_ack, ...
            'Q_ack', Q_ack, ...
            'wn_Error', wn_percentage_error, ...
            'Q_Error', Q_percentage_error, ...
            'R', R, ...
            'C1', C1, ...
            'C2', C2, ...
            'R_rounded', R_rounded, ...
            'C1_rounded', C1_rounded, ...
            'C2_rounded', C2_rounded, ...
            'R_Error', R_pecentage_error, ...
            'C1_Error', C1_pecentage_error, ...
            'C2_Error', C2_pecentage_error, ...
            'Avg_Error', Average_error...
        )];
    end
    
    % This generates strings for LTSpice to look at the response of each
    % set of components.
    % Input '.step param index 1 24 1' (.step param <name> <start> <stop>
    % <increment>) as a spice directive, followed by the
    % generated strings. It creates one string per component.
    r1_table = ".param R1 = table(index";
    r2_table = ".param R2 = table(index";
    r3_table = ".param R3 = table(index";
    c1_table = ".param C1 = table(index";
    c2_table = ".param C2 = table(index";
    
    for i = 1:length(results)
        r1 = results(i).R_rounded;
        r2 = results(i).R_rounded;
        r3 = results(i).R_rounded;
        c1 = results(i).C1_rounded;
        c2 = results(i).C2_rounded;
    
        r1_table = sprintf('%s,%d,%.3g', r1_table, i-1, r1);
        r2_table = sprintf('%s,%d,%.3g', r2_table, i-1, r2);
        r3_table = sprintf('%s,%d,%.3g', r3_table, i-1, r3);
        c1_table = sprintf('%s,%d,%.3g', c1_table, i-1, c1);
        c2_table = sprintf('%s,%d,%.3g', c2_table, i-1, c2);
    end
    
    % Print for copy-paste into LTSpice
    if print_ltspice == true
        fprintf('%s)\n', r1_table);
        fprintf('%s)\n', r2_table);
        fprintf('%s)\n', r3_table);
        fprintf('%s)\n', c1_table);
        fprintf('%s)\n\n', c2_table);
    end
    
    % Display the results as a table, useful to decide what components to use
    if print_table == true
        disp(struct2table(results))
    end
end

function rounded = round_to_E12(value)
% Finds the nearest E12 value for the input component value.
% Generates all the E12 values from pico to Mega.
% Finds the index of the generated values with the smallest difference to
% the input value, then uses the index to return the rounded value.
    E12 = [1.0 1.2 1.5 1.8 2.2 2.7 3.3 3.9 4.7 5.6 6.8 8.2];
    decades = -12:6;
    values = [];
    for d = decades
        values = [values, E12 * 10^d];
    end
    [~, index] = min(abs(values - value));
    rounded = values(index);
end

function err = percent_error(actual, approx)
    err = (actual - approx) / approx * 100;
    err = abs(err);
end
