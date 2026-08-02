function results_table = calliper_csv_results_table(filename, folder)
    % ── READ DATA ───────────────────────────────────────────────
    filepath = fullfile(folder, filename);
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
    
    % ── PARSE METADATA ─────────────────────────────────────────
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
end
