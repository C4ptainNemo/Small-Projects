%% ============================================================
%  Reads the result summary CSV and plots every data column
%  against Piston Force as configurable subplots.
%  Metadata (excluding DP and Piston Force) is shown as an
%  annotation box on the figure.
clear; close all; clc;
%% ============================================================

% ── FILE ─────────────────────────────────────────────────────
filename = 'Result_Summary_PistonForce_3000N_M6Bolt_0N_M4Bolt_0N_2026-06-24_21-35-55.csv';
folder = 'csv_data';
filepath = fullfile(folder, filename);

% ── FONT CONTROLS ────────────────────────────────────────────
fontName        = 'Arial';      % Font family for all text
titleFontSize   = 9;            % Subplot title font size (pt)
axisLabelSize   = 8;            % X/Y axis label font size (pt)
tickLabelSize   = 7;            % Axis tick label font size (pt)
metaFontSize    = 8;            % Metadata annotation font size (pt)
legendFontSize  = 7;            % Legend font size (pt, if used)

% ── SUBPLOT LAYOUT CONTROLS ──────────────────────────────────
% Outer figure dimensions [width, height] in centimetres
figWidthCm      = 50;
figHeightCm     = 35;

% Gaps between and around subplots (normalised 0-1 units)
leftMargin      = 0.06;   % left edge of subplot grid
rightMargin     = 0.78;   % right edge of subplot grid (leave room for metadata)
bottomMargin    = 0.07;   % bottom edge of subplot grid
topMargin       = 0.93;   % top edge of subplot grid
hGap            = 0.04;   % horizontal gap between columns
vGap            = 0.04;   % vertical gap between rows

% ── LINE STYLE CONTROLS ──────────────────────────────────────
lineWidth       = 1.2;
lineColor       = [0.15 0.35 0.65];   % RGB (0-1)
markerStyle     = 'none';             % 'o', 's', 'none', etc.
markerSize      = 4;

% ── METADATA BOX POSITION (normalised figure units) ──────────
metaBoxLeft     = 0.80;   % left edge of annotation box
metaBoxBottom   = 0.10;   % bottom edge
metaBoxWidth    = 0.18;   % width
metaBoxHeight   = 0.80;   % height

%% ── READ DATA ───────────────────────────────────────────────
fid = fopen(filepath, 'r');
if fid == -1
    error('Cannot open file: %s', filepath);
end

% Row 1 – column names
headerLine   = fgetl(fid);
colNames     = strsplit(headerLine, ',');

% Row 2 – quantity labels (e.g. "Maximum", "Working Load")
quantityLine = fgetl(fid);
quantities   = strsplit(quantityLine, ',');

% Row 3 – units
unitLine     = fgetl(fid);
units        = strsplit(unitLine, ',');

% Remaining lines: numeric data until a blank line, then metadata
dataLines = {};
metaLines = {};
inMeta    = false;

while ~feof(fid)
    line = fgetl(fid);
    if ischar(line)
        trimmed = strtrim(line);
        if isempty(trimmed)
            inMeta = true;   % blank line separates data from metadata
            continue;
        end
        if inMeta
            metaLines{end+1} = trimmed; %#ok<AGROW>
        else
            dataLines{end+1} = trimmed; %#ok<AGROW>
        end
    end
end
fclose(fid);

% Parse numeric data matrix
nRows = numel(dataLines);
nCols = numel(colNames);
data  = zeros(nRows, nCols);
for r = 1:nRows
    vals       = strsplit(dataLines{r}, ',');
    data(r, :) = str2double(vals);
end

% X axis: first column (Piston Force)
xData     = data(:, 1);
xLabel    = sprintf('%s (%s)', colNames{1}, units{1});

% Y columns: everything after the first
yData     = data(:, 2:end);
yColNames = colNames(2:end);
yUnits    = units(2:end);
yQty      = quantities(2:end);

%% ── PARSE METADATA ─────────────────────────────────────────
% Exclude lines whose key starts with "DP" or "Piston Force"
excludeKeys = {'DP', 'Piston Force'};
metaDisplay = {};

for i = 1:numel(metaLines)
    parts = strsplit(metaLines{i}, ',');
    if numel(parts) >= 2
        key = strtrim(parts{1});
        val = strtrim(parts{2});
        skip = false;
        for e = 1:numel(excludeKeys)
            if strncmpi(key, excludeKeys{e}, numel(excludeKeys{e}))
                skip = true;
                break;
            end
        end
        if ~skip
            metaDisplay{end+1} = sprintf('%s: %s', key, val); %#ok<AGROW>
        end
    end
end

%% ── FIGURE & LAYOUT ────────────────────────────────────────
nPlots = size(yData, 2);
nCols_ = ceil(sqrt(nPlots));
nRows_ = ceil(nPlots / nCols_);

% Convert cm to inches for MATLAB figure size
fig = figure('Units', 'centimeters', ...
             'Position', [2 2 figWidthCm figHeightCm], ...
             'Color', 'white');

% Compute individual subplot size in normalised units
totalW    = rightMargin - leftMargin;
totalH    = topMargin   - bottomMargin;
subW      = (totalW - hGap * (nCols_ - 1)) / nCols_;
subH      = (totalH - vGap * (nRows_ - 1)) / nRows_;

axHandles = gobjects(nPlots, 1);

for p = 1:nPlots
    row = ceil(p / nCols_);   % 1 = top
    col = mod(p - 1, nCols_) + 1;

    % Compute position [left bottom width height] in normalised units
    posL = leftMargin  + (col - 1) * (subW + hGap);
    posB = topMargin   - row * subH - (row - 1) * vGap;   % from top
    pos  = [posL, posB, subW, subH];

    ax = axes('Parent', fig, 'Position', pos); %#ok<LAXES>
    axHandles(p) = ax;

    % Build Y-axis label: "Quantity (unit)"
    yLbl = sprintf('%s (%s)', yQty{p}, yUnits{p});

    plot(ax, xData, yData(:, p), ...
        'Color',      lineColor, ...
        'LineWidth',  lineWidth, ...
        'Marker',     markerStyle, ...
        'MarkerSize', markerSize);

    % Title uses the column name (remove leading group prefix for brevity)
    title(ax, yColNames{p}, ...
        'FontName', fontName, ...
        'FontSize', titleFontSize, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none');

    xlabel(ax, xLabel, ...
        'FontName', fontName, ...
        'FontSize', axisLabelSize);

    ylabel(ax, yLbl, ...
        'FontName', fontName, ...
        'FontSize', axisLabelSize);

    set(ax, ...
        'FontName',  fontName, ...
        'FontSize',  tickLabelSize, ...
        'Box',       'on', ...
        'XGrid',     'on', ...
        'YGrid',     'on', ...
        'GridAlpha', 0.25);
end

% Hide any unused axes (if nPlots not a perfect grid)
for p = nPlots + 1 : nRows_ * nCols_
    row = ceil(p / nCols_);
    col = mod(p - 1, nCols_) + 1;
    posL = leftMargin  + (col - 1) * (subW + hGap);
    posB = topMargin   - row * subH - (row - 1) * vGap;
    ax   = axes('Parent', fig, 'Position', [posL, posB, subW, subH], ...
                'Visible', 'off'); %#ok<LAXES>
end

%% ── METADATA ANNOTATION ────────────────────────────────────
% Draw a box and place text inside it
annotation(fig, 'rectangle', ...
    [metaBoxLeft, metaBoxBottom, metaBoxWidth, metaBoxHeight], ...
    'Color', [0.5 0.5 0.5], ...
    'LineWidth', 0.8, ...
    'FaceColor', [0.97 0.97 0.97]);

metaStr = sprintf('\\bfSimulation Metadata\\rm\n\n%s', ...
                  strjoin(metaDisplay, '\n'));

annotation(fig, 'textbox', ...
    [metaBoxLeft + 0.005, metaBoxBottom + 0.01, ...
     metaBoxWidth - 0.01, metaBoxHeight - 0.02], ...
    'String',          metaStr, ...
    'FontName',        fontName, ...
    'FontSize',        metaFontSize, ...
    'EdgeColor',       'none', ...
    'BackgroundColor', 'none', ...
    'VerticalAlignment',   'top', ...
    'HorizontalAlignment', 'left', ...
    'Interpreter',     'tex', ...
    'FitBoxToText',    'off');

%% ── SAVE FIGURE ───────────────────────────────────
% Uncomment to save as PNG or PDF:
[folder, name, ~] = fileparts(filepath);          % folder and name from filepath
figure_filename = sprintf('%s_plots', name);

% Ensure output full paths in the same folder as the data file
pngPath = fullfile(folder, [figure_filename, '.png']);
pdfPath = fullfile(folder, [figure_filename, '.pdf']);

% Save PNG and PDF using the derived filename in the same folder
exportgraphics(fig, pngPath, 'Resolution', 300);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');