function [modifiedStruct, nanStatsStruct] = modifyStructWithGap7_sx3(input)
% MODIFYSTRUCTWITHGAP7_SX3 Filters outliers and handles temporal gaps with flexible alignment
%
% INPUT:
%   input struct with fields:
%       - sogliaiqr_up      : upper IQR multiplier
%       - sogliaiqr_do      : lower IQR multiplier
%       - originalStruct    : struct with tables (first column: datetime, others: values)
%       - gap               : step in days between generated dates
%       - flexibleAlignment : boolean, true for flexible alignment
%       - n                 : day tolerance for flexible alignment
%
% OUTPUT:
%   - modifiedStruct : struct with filtered tables and missing-value placeholders
%   - nanStatsStruct : struct with Before/After statistics + outlier counts

    sogliaiqr_up = input.sogliaiqr_up;
    sogliaiqr_do = input.sogliaiqr_do;
    originalStruct = input.originalStruct;
    gap = input.gap;
    flexibleAlignment = isfield(input, 'flexibleAlignment') && input.flexibleAlignment;
    n = 1;
    if flexibleAlignment && isfield(input, 'n')
        n = input.n;
    end

    modifiedStruct = struct();
    nanStatsStruct = struct();
    fields = fieldnames(originalStruct);

    for i = 1:length(fields)
        fieldName = fields{i};
        originalTable = originalStruct.(fieldName);
        originalDates = originalTable{:, 1};
        originalValues = table2array(originalTable(:, 2:end));
        colNames = originalTable.Properties.VariableNames(2:end);

        % Find valid dates
        validRows = any(~isnan(originalValues), 2);
        validDates = originalDates(validRows);
        if isempty(validDates)
            warning(['No valid dates found in ' fieldName ', all NaN.']);
            modifiedStruct.(fieldName) = originalTable;
            nanStatsStruct.(fieldName) = [];
            continue;
        end

        % ==========================================================
        % 1) OUTLIER REMOVAL
        % ==========================================================
        filteredValues = originalValues;
        outlierCount = zeros(1, size(filteredValues, 2));
        for col = 1:size(filteredValues, 2)
            colData = filteredValues(:, col);
            med = median(colData, 'omitnan');
            iqrVal = iqr(colData);
            up = med + sogliaiqr_up * iqrVal;
            low = med - sogliaiqr_do * iqrVal;
            outliers = colData > up | colData < low;
            outlierCount(col) = sum(outliers);
            filteredValues(outliers, col) = NaN;
        end

        % ==========================================================
        % 2) Build filtered table after outlier removal
        % ==========================================================
        filteredTable = [table(originalDates, 'VariableNames', {'Date'}), ...
                         array2table(filteredValues, 'VariableNames', colNames)];

        % ==========================================================
        % 3) BEFORE STATISTICS (after outlier removal)
        % ==========================================================
        totalBefore = size(filteredValues, 1);
        nonNaNBefore = sum(~isnan(filteredValues), 1);
        nanBefore = totalBefore - nonNaNBefore;
        percentNanBefore = 100 * nanBefore / totalBefore;
        percentNonNanBefore = 100 * nonNaNBefore / totalBefore;

        maxBefore = max(filteredValues, [], 1, 'omitnan');
        minBefore = min(filteredValues, [], 1, 'omitnan');

        statsBefore = table(nonNaNBefore', nanBefore', percentNonNanBefore', percentNanBefore', ...
                            maxBefore', minBefore', outlierCount', ...
            'VariableNames', {'NonNaN', 'NaN', 'PercentNonNaN', 'PercentNaN', ...
                              'Max', 'Min', 'OutlierRemoved'}, ...
            'RowNames', colNames);

        % ==========================================================
        % 4) APPLY ALIGNMENT
        % ==========================================================
        [newDates, newValues] = divide_myTableFlexible_sx(filteredTable, gap, n);

        % ==========================================================
        % 5) AFTER STATISTICS
        % ==========================================================
        totalAfter = size(newValues, 1);
        nonNaNAfter = sum(~isnan(newValues), 1);
        nanAfter = totalAfter - nonNaNAfter;
        percentNanAfter = 100 * nanAfter / totalAfter;
        percentNonNanAfter = 100 * nonNaNAfter / totalAfter;

        maxAfter = max(newValues, [], 1, 'omitnan');
        minAfter = min(newValues, [], 1, 'omitnan');

        statsAfter = table(nonNaNAfter', nanAfter', percentNonNanAfter', percentNanAfter', ...
                           maxAfter', minAfter', ...
            'VariableNames', {'NonNaN', 'NaN', 'PercentNonNaN', 'PercentNaN', ...
                              'Max', 'Min'}, ...
            'RowNames', colNames);

        % ==========================================================
        % 6) SAVE STATISTICS
        % ==========================================================
        nanStatsStruct.(fieldName).Before = statsBefore;
        nanStatsStruct.(fieldName).After = statsAfter;
        nanStatsStruct.(fieldName).Outliers = outlierCount;

        % ==========================================================
        % 7) FINAL TABLE
        % ==========================================================
        modifiedTable = [table(newDates, 'VariableNames', {'Date'}), ...
                         array2table(newValues, 'VariableNames', colNames)];
        modifiedStruct.(fieldName) = modifiedTable;
    end

    disp('Processing complete: outliers removed, multi-level statistics saved.');
end


%% ════════════════════════════════════════════════════════════════════════
%%  LOCAL FUNCTION
%% ════════════════════════════════════════════════════════════════════════

function [newDates, newVals] = divide_myTableFlexible_sx(inputTable, gap, n)
% DIVIDE_MYTABLEFLEXIBLE_SX  Resamples a time table to a regular interval (gap days).
% Left-side alignment: each value is assigned to the nearest target date
% on the left (i.e. preceding or equal) within a tolerance of n days.
%
% INPUT:
%   - inputTable : table with first column datetime, remaining columns numeric
%   - gap        : step size in days between generated dates
%   - n          : day tolerance for matching (looks back up to n days)
%
% OUTPUT:
%   - newDates   : vector of regular target dates
%   - newVals    : matrix of aligned values (NaN where no match found)

    dates  = inputTable{:, 1};
    values = table2array(inputTable(:, 2:end));

    % Keep only rows with at least one valid value
    maskValid  = any(~isnan(values), 2);
    validDates = dates(maskValid);

    if isempty(validDates)
        newDates = datetime.empty;
        newVals  = [];
        return;
    end

    startDate = min(validDates);
    endDate   = max(validDates);

    newDates = (startDate : days(gap) : endDate)';
    if newDates(end) < endDate
        newDates = [newDates; endDate];
    end

    nDates = length(newDates);
    nVars  = size(values, 2);
    newVals = NaN(nDates, nVars);

    % Left-side alignment: find the closest observation <= target within n days
    for i = 1:nDates
        dateTarget = newDates(i);

        % Signed differences: positive = future, negative = past
        diffs = days(dates - dateTarget);

        % Candidates: dates on or before target, no more than n days back
        candidateIdx = find(diffs <= 0 & diffs >= -n);

        if ~isempty(candidateIdx)
            % Take the last candidate (closest to the left of target)
            idx = candidateIdx(end);
            newVals(i, :) = values(idx, :);
        end
    end
end
