%% ============================================================
%  Reads multiple result summary CSVs and plots every data column
%  against Piston Force as configurable subplots.
%  All files are overlaid on the same figures, distinguished by
%  colour and line style. Metadata and legend occupy the final
%  two subplot cells.
clear; close all; clc;
%% ============================================================

% ── FILES ────────────────────────────────────────────────────
folder = 'csv_data';
filenames = {
    'Result_Summary_PistonForce_9N_M6Bolt_0N_M4Bolt_0N.csv';
    'Result_Summary_PistonForce_9N_M6Bolt_1000N_M4Bolt_500N.csv';
    'Result_Summary_PistonForce_9N_M6Bolt_3000N_M4Bolt_1500N.csv';
    'Result_Summary_PistonForce_9N_M6Bolt_5000N_M4Bolt_2500N.csv';
    'Result_Summary_PistonForce_9N_M6Bolt_7000N_M4Bolt_3500N.csv';
    'Result_Summary_PistonForce_9N_M6Bolt_10000N_M4Bolt_5000N.csv';
};

% ── LEGEND LABELS (one per file) ─────────────────────────────
% Ensure these are in the same order as the filenames
legendLabels = {
    '0 / 0';
    '1000 / 500';
    '3000 / 1500'
    '5000 / 2500'
    '7000 / 3500'
    '10000 / 5000';
};

% ── FONT CONTROLS ────────────────────────────────────────────
fontName        = 'Arial';
titleFontSize   = 9;
axisLabelSize   = 8;
tickLabelSize   = 7;
metaFontSize    = 8;
legendFontSize  = 8;

% ── SUBPLOT LAYOUT CONTROLS ──────────────────────────────────
figWidthCm      = 50;
figHeightCm     = 35;

leftMargin      = 0.06;
rightMargin     = 0.97;   % now uses full width
bottomMargin    = 0.07;
topMargin       = 0.93;
hGap            = 0.04;
vGap            = 0.04;

% ── LINE STYLE CONTROLS ──────────────────────────────────────
lineWidth       = 1.2;
markerSize      = 3;

lineColors = [
    0.15 0.35 0.65;
    0.85 0.20 0.10;
    0.10 0.60 0.20;
    0.80 0.50 0.00;
    0.50 0.10 0.70;
    0.00 0.60 0.70;
];
lineStyles   = {'-', '--', ':', '-.'};
markerStyles = {'none', 'o', 's', '^', 'd', 'v'};

%% ── LOAD ALL FILES ──────────────────────────────────────────
nFiles    = numel(filenames);
allData   = cell(nFiles, 1);
allMeta   = cell(nFiles, 1);
colNames  = {};
quantities = {};
units     = {};

for f = 1:nFiles
    filepath = fullfile(folder, filenames{f});
    fid = fopen(filepath, 'r');
    if fid == -1
        error('Cannot open file: %s', filepath);
    end

    headerLine   = fgetl(fid);
    quantityLine = fgetl(fid);
    unitLine     = fgetl(fid);

    if f == 1
        colNames   = strsplit(headerLine,   ',');
        quantities = strsplit(quantityLine, ',');
        units      = strsplit(unitLine,     ',');
    end

    dataLines = {};
    metaLines = {};
    inMeta    = false;

    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line)
            trimmed = strtrim(line);
            if isempty(trimmed)
                inMeta = true;
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

    nRows   = numel(dataLines);
    nCols   = numel(colNames);
    dataMat = zeros(nRows, nCols);
    for r = 1:nRows
        vals          = strsplit(dataLines{r}, ',');
        dataMat(r, :) = str2double(vals);
    end

    allData{f} = dataMat;
    allMeta{f} = metaLines;
end

%% ── PARSE METADATA (first file only, Solution Time listed per file) ───
excludeKeys    = {'DP', 'Piston Force', 'M6 Pretension', 'M4 Pretension', 'Solution Date'};
solutionTimes  = cell(nFiles, 1);   % collect Solution Time per file
metaDisplay    = {};                 % shared metadata from file 1

for f = 1:nFiles
    metaLinesF = allMeta{f};
    for i = 1:numel(metaLinesF)
        parts = strsplit(metaLinesF{i}, ',');
        if numel(parts) >= 2
            key = strtrim(parts{1});
            val = strtrim(parts{2});

            % Always capture Solution Time for every file
            if strcmpi(key, 'Solution Time')
                solutionTimes{f} = val;
                continue;
            end

            % For file 1 only, collect remaining metadata (excluding DP / Piston Force)
            if f == 1
                skip = false;
                for e = 1:numel(excludeKeys)
                    if strncmpi(key, excludeKeys{e}, numel(excludeKeys{e}))
                        skip = true; break;
                    end
                end
                if ~skip
                    metaDisplay{end+1} = sprintf('%s: %s', key, val); %#ok<AGROW>
                end
            end
        end
    end
end

%% ── LEGEND LABELS ───────────────────────────────────────────
if isempty(legendLabels)
    legendLabels = cell(nFiles, 1);
    for f = 1:nFiles
        [~, name, ~]    = fileparts(filenames{f});
        legendLabels{f} = strrep(name, '_', ' ');
    end
end

%% ── FIGURE & LAYOUT ─────────────────────────────────────────
xLabel    = sprintf('%s (%s)', colNames{1}, units{1});
yColNames = colNames(2:end);
yUnits    = units(2:end);
yQty      = quantities(2:end);
nPlots    = numel(yColNames);

% Total cells = data plots + legend cell + metadata cell
nCells = nPlots + 2;
nCols_ = ceil(sqrt(nCells));
nRows_ = ceil(nCells / nCols_);

fig = figure('Units', 'centimeters', ...
             'Position', [2 2 figWidthCm figHeightCm], ...
             'Color', 'white');

totalW = rightMargin - leftMargin;
totalH = topMargin   - bottomMargin;
subW   = (totalW - hGap * (nCols_ - 1)) / nCols_;
subH   = (totalH - vGap * (nRows_ - 1)) / nRows_;

lineHandles = gobjects(nFiles, 1);

% Helper: get [left bottom width height] for cell index p
getCellPos = @(p) [...
    leftMargin + (mod(p-1, nCols_))       * (subW + hGap), ...
    topMargin  - ceil(p / nCols_) * subH  - (ceil(p / nCols_) - 1) * vGap, ...
    subW, subH];

%% ── DATA SUBPLOTS ───────────────────────────────────────────
for p = 1:nPlots
    ax = axes('Parent', fig, 'Position', getCellPos(p)); %#ok<LAXES>
    hold(ax, 'on');

    yLbl = sprintf('%s (%s)', yQty{p}, yUnits{p});

    for f = 1:nFiles
        cIdx = mod(f-1, size(lineColors,  1)) + 1;
        lIdx = mod(f-1, numel(lineStyles))    + 1;
        mIdx = mod(f-1, numel(markerStyles))  + 1;

        xf = allData{f}(:, 1);
        yf = allData{f}(:, p + 1);

        h = plot(ax, xf, yf, ...
            'Color',      lineColors(cIdx, :), ...
            'LineStyle',  lineStyles{lIdx}, ...
            'LineWidth',  lineWidth, ...
            'Marker',     markerStyles{mIdx}, ...
            'MarkerSize', markerSize);

        if p == 1
            lineHandles(f) = h;
        end
    end

    hold(ax, 'off');

    title(ax, yColNames{p}, ...
        'FontName', fontName, 'FontSize', titleFontSize, ...
        'FontWeight', 'bold', 'Interpreter', 'none');
    xlabel(ax, xLabel,  'FontName', fontName, 'FontSize', axisLabelSize);
    ylabel(ax, yLbl,    'FontName', fontName, 'FontSize', axisLabelSize);
    set(ax, 'FontName', fontName, 'FontSize', tickLabelSize, ...
        'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.25);
end

%% ── LEGEND SUBPLOT (second-to-last cell) ────────────────────
axLeg = axes('Parent', fig, 'Position', getCellPos(nPlots + 1)); %#ok<LAXES>
set(axLeg, 'Visible', 'off');

% Dummy plots to attach legend to this axes
hold(axLeg, 'on');
dummyH = gobjects(nFiles, 1);
for f = 1:nFiles
    cIdx = mod(f-1, size(lineColors,  1)) + 1;
    lIdx = mod(f-1, numel(lineStyles))    + 1;
    mIdx = mod(f-1, numel(markerStyles))  + 1;
    dummyH(f) = plot(axLeg, NaN, NaN, ...
        'Color',     lineColors(cIdx, :), ...
        'LineStyle', lineStyles{lIdx}, ...
        'LineWidth', lineWidth, ...
        'Marker',    markerStyles{mIdx}, ...
        'MarkerSize', markerSize);
end
hold(axLeg, 'off');

lgd = legend(axLeg, dummyH, legendLabels, ...
    'FontName',        fontName, ...
    'FontSize',        legendFontSize, ...
    'Box',             'on', ...
    'Location',        'northwest');
title(lgd, '\bfLegend: M6/M4 Pretension (N)', 'FontName', fontName, 'FontSize', legendFontSize, ...
    'Interpreter', 'tex');

%% ── METADATA SUBPLOT (last cell) ────────────────────────────
axMeta = axes('Parent', fig, 'Position', getCellPos(nPlots + 2)); %#ok<LAXES>

set(axMeta, ...
    'XLim',     [0 1], ...
    'YLim',     [0 1], ...
    'XTick',    [], ...
    'YTick',    [], ...
    'Box',      'on', ...
    'LineWidth', 0.8, ...
    'XColor',   [0.5 0.5 0.5], ...
    'YColor',   [0.5 0.5 0.5], ...
    'Color',    [0.97 0.97 0.97]);

% Build Solution Time list string, one entry per file
solTimeLines = cell(nFiles, 1);
for f = 1:nFiles
    [~, fname, ~]     = fileparts(filenames{f});
    solTimeLines{f}   = sprintf('  %s: %s s', legendLabels{f}, solutionTimes{f});
end
solTimeStr = strjoin(solTimeLines, '\n');

metaStr = sprintf('\\bfSimulation Metadata\\rm\n\n%s\n\n\\bfSolution Time\\rm\n%s', ...
                  strjoin(metaDisplay, '\n'), solTimeStr);

text(axMeta, 0.05, 0.95, metaStr, ...
    'Units',               'normalized', ...
    'FontName',            fontName, ...
    'FontSize',            metaFontSize, ...
    'VerticalAlignment',   'top', ...
    'HorizontalAlignment', 'left', ...
    'Interpreter',         'tex');

%% ── SAVE FIGURE ─────────────────────────────────────────────
[~, firstName, ~]   = fileparts(filenames{1});
%figure_filename     = sprintf('%s_plots_combined', firstName);
figure_filename     = sprintf('Combined_Results_Plots');

pngPath = fullfile(folder, [figure_filename, '.png']);
pdfPath = fullfile(folder, [figure_filename, '.pdf']);

exportgraphics(fig, pngPath, 'Resolution', 600);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');