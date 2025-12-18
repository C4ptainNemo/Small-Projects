function results = component_calculator_vcvs ...
                            (Q, wn, C_base, print_table, print_ltspice)
% Function for component tuning for VCVS filter.
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
        % Taken from the lecture notes by Dr Jasmine Banks
        m = 10^(k / 24);
        n = ((Q * (1 + m)) / sqrt(m))^2;
        C1 = C_base;
        C2 = n * C1;
        R1 = 1 / (wn * C1 * sqrt(m * n));
        R2 = m * R1;

        % Round to E12
        R1_rounded = round_to_E12(R1);
        R2_rounded = round_to_E12(R2);
        C1_rounded = round_to_E12(C1);
        C2_rounded = round_to_E12(C2);

        wn_ack = 1 / sqrt(R1_rounded*R2_rounded*C1_rounded*C2_rounded);
        Q_ack = sqrt(R1_rounded*R2_rounded*C1_rounded*C2_rounded) ...
                / (R2_rounded*C1_rounded + R1_rounded*C1_rounded);
        wn_percentage_error = percent_error(wn, wn_ack);
        Q_percentage_error = percent_error(Q, Q_ack);
        avg_wn_Q_error = (wn_percentage_error + Q_percentage_error) / 2;

        R1_pecentage_error = percent_error(R1, R1_rounded);
        R2_pecentage_error = percent_error(R2, R2_rounded);
        C1_pecentage_error = percent_error(C1, C1_rounded);
        C2_pecentage_error = percent_error(C2, C2_rounded);
        Average_error = ...
                (R1_pecentage_error + R2_pecentage_error + C2_pecentage_error) / 3;

        % Store results
        results = [results; struct(...
            'k', k, ...
            'wn_ack', wn_ack, ...
            'Q_ack', Q_ack, ...
            'wn_err', wn_percentage_error, ...
            'Q_err', Q_percentage_error, ...
            'Avg_wn_Q_err', avg_wn_Q_error, ...
            'R1', R1, ...
            'R2', R2, ...
            'C1', C1, ...
            'C2', C2, ...
            'R1_rnd', R1_rounded, ...
            'R2_rnd', R2_rounded, ...
            'C1_rnd', C1_rounded, ...
            'C2_rnd', C2_rounded, ...
            'R1_err', R1_pecentage_error, ...
            'R2_err', R2_pecentage_error, ...
            'C2_err', C2_pecentage_error, ...
            'Avg_err', Average_error...
        )];
    end

    % This generates strings for LTSpice to look at the response for each
    % set of components.
    % Input '.step param index 1 24 1' (.step param <name> <start> <stop>
    % <increment>) as a spice directive, followed by the
    % generated strings. It creates one string per component.
    r1_table = ".param R1 = table(index";
    r2_table = ".param R2 = table(index";
    c1_table = ".param C1 = table(index";
    c2_table = ".param C2 = table(index";
    for i = 1:length(results)
        r1 = results(i).R1_rnd;
        r2 = results(i).R2_rnd;
        c1 = results(i).C1_rnd;
        c2 = results(i).C2_rnd;
    
        r1_table = sprintf('%s,%d,%.3g', r1_table, i-1, r1);
        r2_table = sprintf('%s,%d,%.3g', r2_table, i-1, r2);
        c1_table = sprintf('%s,%d,%.3g', c1_table, i-1, c1);
        c2_table = sprintf('%s,%d,%.3g', c2_table, i-1, c2);
    end
    
    % Print for copy-paste into LTSpice
    if print_ltspice == true
        fprintf('%s)\n', r1_table);
        fprintf('%s)\n', r2_table);
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
% Finds the index of the generated value with the smallest difference to
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
    err = abs(actual - approx) / actual * 100;
end