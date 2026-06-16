function [H, C] = complexity_entropy_universal(data, dx, dy, taux, tauy, normalize, use_miller,probs, tie_precision)
% COMPLEXITY_ENTROPY_NAN
% Entropia di permutazione H e complessità C per dati 1D/2D con NaN
% Con opzioni: normalizzazione e correzione Miller-Madow attivabile

    if nargin < 2, dx = 3; end
    if nargin < 3, dy = 1; end
    if nargin < 4, taux = 1; end
    if nargin < 5, tauy = 1; end
    if nargin < 6 || isempty(normalize), normalize = true; end
    if nargin < 7 || isempty(use_miller), use_miller = false; end
    if nargin < 8, probs = false; end
    if nargin < 9, tie_precision = []; end

    n = factorial(dx*dy);

    % --- Distribuzione ---
    if ~probs
        [~, probabilities, N_valid] = ordinal_distribution_nan(data, dx, dy, taux, tauy, false, tie_precision, false);
    else
        probabilities = data(:);
        probabilities = probabilities(probabilities > 0);
        N_valid = [];
    end

    % --- Entropia ---
    H = permutation_entropy_nan(probabilities, dx, dy, taux, tauy, 2, normalize, true, tie_precision, N_valid, use_miller);

    % --- Distribuzione uniforme ---
    U = ones(n,1)/n;
    P_full = zeros(n,1);
    P_full(1:numel(probabilities)) = probabilities;

    % --- Jensen-Shannon ---
    P_plus_U_over_2 = (P_full + U)/2;
    S_pu = -sum(P_plus_U_over_2 .* log(P_plus_U_over_2 + eps));
    S_p  = -0.5 * sum(P_full .* log(P_full + eps));
    S_u  = -0.5 * sum(U .* log(U + eps));
    js_div = S_pu - S_p - S_u;

    % --- JS max ---
    js_div_max = -0.5 * ( ((n+1)/n)*log(n+1) + log(n) - 2*log(2*n) );

    % --- Complessità ---
    C = H * js_div / js_div_max;
end

%% ----------------------------------------------------------------------
function H = permutation_entropy_nan(data, dx, dy, taux, tauy, base, normalized, probs, tie_precision, N_valid, use_miller)
% Entropia di permutazione con Miller opzionale

    if nargin < 2, dx = 3; end
    if nargin < 3, dy = 1; end
    if nargin < 4, taux = 1; end
    if nargin < 5, tauy = 1; end
    if nargin < 6, base = 2; end
    if nargin < 7 || isempty(normalized), normalized = true; end
    if nargin < 8, probs = false; end
    if nargin < 9, tie_precision = []; end
    if nargin < 10, N_valid = []; end
    if nargin < 11 || isempty(use_miller), use_miller = true; end

    if ~probs
        [~, probabilities, N_auto] = ordinal_distribution_nan(data, dx, dy, taux, tauy, false, tie_precision, false);
        if isempty(N_valid), N_valid = N_auto; end
    else
        probabilities = data(:);
        probabilities = probabilities(probabilities > 0);
    end

    % log base
    if isequal(base,'e'), logfunc = @(x) log(x);
    elseif isequal(base,2), logfunc = @(x) log2(x);
    else, error('Base must be ''e'' or 2'); end

    % Entropia osservata
    H_obs = -sum(probabilities .* logfunc(probabilities));

    % --- Miller-Madow opzionale ---
    bias = 0;
    if use_miller && ~isempty(N_valid) && N_valid > 0
        K = sum(probabilities > 0);
        if isequal(base,2)
            bias = (K - 1) / (2 * N_valid * log(2));
        else
            bias = (K - 1) / (2 * N_valid);
        end
    end

    H_abs = H_obs + bias;

    % --- Normalizzazione opzionale ---
    if normalized
        H_max = logfunc(factorial(dx*dy));
        H = H_abs / H_max;
    else
        H = H_abs;
    end
end

%% ----------------------------------------------------------------------
function [symbols, probabilities, N_valid] = ordinal_distribution_nan(data, dx, dy, taux, tauy, return_missing, tie_precision, ordered)

    if nargin < 2, dx = 3; end
    if nargin < 3, dy = 1; end
    if nargin < 4, taux = 1; end
    if nargin < 5, tauy = 1; end
    if nargin < 6, return_missing = false; end
    if nargin < 7, tie_precision = []; end
    if nargin < 8, ordered = false; end

    symbols = ordinal_sequence_nan(data, dx, dy, taux, tauy, true, tie_precision);
    if ndims(symbols) == 3, symbols = reshape(symbols, [], dx*dy); end

    valid_idx = all(~isnan(symbols), 2);
    symbols = symbols(valid_idx,:);
    N_valid = size(symbols,1);

    if N_valid == 0
        probabilities = [];
        symbols = [];
        return;
    end

    [unique_symbols, ~, ic] = unique(symbols,'rows');
    counts = accumarray(ic,1);
    probabilities = counts / N_valid;

    if return_missing
        total_patterns = factorial(dx*dy);
        if size(unique_symbols,1) < total_patterns
            all_symbols = flipud(perms(0:(dx*dy-1)));
            [~, ia] = setdiff(all_symbols, unique_symbols,'rows','stable');
            missing_symbols = all_symbols(ia,:);
            symbols = [unique_symbols; missing_symbols];
            probabilities = [probabilities; zeros(size(missing_symbols,1),1)];
        else
            symbols = unique_symbols;
        end
    else
        symbols = unique_symbols;
    end
end

%% ----------------------------------------------------------------------
function symbols = ordinal_sequence_nan(data, dx, dy, taux, tauy, overlapping, tie_precision)

    if nargin < 2, dx = 3; end
    if nargin < 3, dy = 1; end
    if nargin < 4, taux = 1; end
    if nargin < 5, tauy = 1; end
    if nargin < 6, overlapping = true; end
    if nargin < 7, tie_precision = []; end

    data = double(data);
    [ny, nx] = size(data);
    if isvector(data), data = data(:)'; ny = 1; nx = length(data); end
    if ~isempty(tie_precision), data = round(data,tie_precision); end

    if ny == 1
        indices = (1:taux:(1+(dx-1)*taux))' + (0:(nx-(dx-1)*taux)-1);
        indices(:,any(indices>nx,1)) = [];
        partitions = data(indices)';
        [~, symbols] = sort(partitions,2,'ascend');
        symbols(any(isnan(partitions),2),:) = NaN;
        symbols = symbols - 1;
    else
        symbols = zeros(ny,nx,dx*dy);
        for i = 1:ny-(dy-1)*tauy
            for j = 1:nx-(dx-1)*taux
                block = zeros(dy,dx);
                for ii=0:(dy-1)
                    for jj=0:(dx-1)
                        block(ii+1,jj+1) = data(i+ii*tauy,j+jj*taux);
                    end
                end
                if any(isnan(block(:))), symbols(i,j,:) = NaN;
                else [~, symbols(i,j,:)] = sort(block(:)); end
            end
        end
        symbols = symbols - 1;
    end
end
