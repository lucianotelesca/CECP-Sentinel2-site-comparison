function result = minimum_complexity_entropy_0(dx, dy, size, plotFlag)
% MINIMUM_COMPLEXITY_ENTROPY
% Calcola la curva teorica di minima complessità statistica
% per una distribuzione di permutazioni dx x dy
% INPUT:
%   dx, dy : dimensione dell'embedding
%   size   : numero di punti sulla curva
%   plotFlag : se true, plotta la curva (default true)
%
% OUTPUT:
%   result : matrice [H, C] della curva minima

size = size + 1;
N = factorial(dx * dy);
prob_params = linspace(1/N, 1, size-1);
uniform_dist = ones(1, N) / N;

hc_ = [];
for i = 1:(size-1)
    probabilities = ones(1, N) * (1 - prob_params(i)) / (N - 1);
    probabilities(1) = prob_params(i);

    % Calcolo entropia
    h = permutation_entropy(probabilities, dx, dy, 1, 1, 2, true, true);

    % Calcolo complessità
    p_plus_u_over_2 = (uniform_dist + probabilities) / 2;
    s_of_p_plus_u_over_2 = -sum(p_plus_u_over_2 .* log(p_plus_u_over_2));

    probabilities_non_zero = probabilities(probabilities ~= 0);
    s_of_p_over_2 = -sum(probabilities_non_zero .* log(probabilities_non_zero)) / 2;
    s_of_u_over_2 = log(N) / 2;

    js_div_max = -0.5 * (((N + 1) / N) * log(N + 1) + log(N) - 2 * log(2 * N));
    js_div = s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2;

    hc_ = [hc_; h, h * js_div / js_div_max];
end

result = flipud(hc_);

% --- Plot della curva minima ---
if plotFlag
    figure; hold on;
    plot(result(:,1), result(:,2), 'r-', 'LineWidth', 2);
    xlabel('Entropia H'); ylabel('Complessità C');
    title(sprintf('Curva minima H-C per dx=%d, dy=%d', dx, dy));
    grid on; box on;
end

end
