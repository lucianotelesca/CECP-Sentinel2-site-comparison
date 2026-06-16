function result = maximum_complexity_entropy_0(dx, dy, m, plotFlag)
% MAXIMUM_COMPLEXITY_ENTROPY
% Calcola la curva teorica di massima complessità statistica
% per una distribuzione di permutazioni dx x dy
%
% INPUT:
%   dx, dy : dimensione dell'embedding
%   m      : numero di punti sulla curva (default 100)
%   plotFlag : se true, plotta la curva (default true)
%
% OUTPUT:
%   result : matrice [H, C] della curva massima teorica

    N = factorial(dx * dy);
    hlist_ = zeros(N - 1, m);
    clist_ = zeros(N - 1, m);

    uniform_dist = ones(1, N) / N;
    prob_params = linspace(0, 1 / N, m);

    for i = 1:N - 1
        p = zeros(1, N);

        for k = 1:length(prob_params)
            p(1) = prob_params(k);
            for j = 2:(N - i)
                p(j) = (1 - prob_params(k)) / (N - i - 1);
            end

            % --- Calcolo entropia ---
            h = permutation_entropy(p, dx, dy, 1, 1, 2, true, true);

            % --- Calcolo complessità ---
            p_plus_u_over_2 = (uniform_dist + p) / 2;
            s_of_p_plus_u_over_2 = -sum(p_plus_u_over_2 .* log(p_plus_u_over_2));

            p_non_zero = p(p > 0);
            s_of_p_over_2 = -sum(p_non_zero .* log(p_non_zero)) / 2;
            s_of_u_over_2 = log(N) / 2;

            js_div_max = -0.5 * (((N + 1) / N) * log(N + 1) + log(N) - 2 * log(2 * N));
            js_div = s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2;

            hlist_(i, k) = h;
            clist_(i, k) = h * js_div / js_div_max;
        end
    end

    % --- Ordinamento e output ---
    hlist_ = hlist_(:);
    clist_ = clist_(:);
    [~, args] = sort(hlist_);
    result = [hlist_(args), clist_(args)];

    % --- Plot della curva massima ---
    if plotFlag
        figure; hold on;
        plot(result(:,1), result(:,2), 'k--', 'LineWidth', 2);
        xlabel('Entropia H'); ylabel('Complessità C');
        title(sprintf('Curva massima H–C per dx=%d, dy=%d', dx, dy));
        grid on; box on;
    end

end
