function results = CECP_compare_sites_allbands(struct_CP, struct_FOL, dx, options)
% CECP_COMPARE_SITES_ALLBANDS
%   Computes H and C for every pixel of two sites (CP and FOL) across all bands,
%   and produces three figures with subplots: CECP, ROC, statistics.
%
% SYNTAX
%   results = CECP_compare_sites_allbands(struct_CP, struct_FOL)
%   results = CECP_compare_sites_allbands(struct_CP, struct_FOL, dx)
%   results = CECP_compare_sites_allbands(struct_CP, struct_FOL, dx, options)
%
% INPUT
%   struct_CP  : infected-site structure  (fields = bands)
%   struct_FOL : healthy-site structure   (fields = bands)
%   dx         : permutation order        (default: 5)
%
% OPTIONS (name-value)
%   dy          : y-pattern size          (default: 1)
%   taux        : x lag                   (default: 1)
%   tauy        : y lag                   (default: 1)
%   normalize   : normalise H             (default: true)
%   use_miller  : Miller-Madow correction (default: false)
%   alpha       : significance level      (default: 0.05)
%   show_curves : CECP limit curves       (default: false)
%   bands       : cell array of bands     (default: all)
%   save_fig    : save figures            (default: true)
%   out_dir     : output folder           (default: pwd)

    arguments
        struct_CP
        struct_FOL
        dx          (1,1) double = 5
        options.dy          (1,1) double  = 1
        options.taux        (1,1) double  = 1
        options.tauy        (1,1) double  = 1
        options.normalize   (1,1) logical = true
        options.use_miller  (1,1) logical = false
        options.show_curves (1,1) logical = false
        options.bands       (1,:) cell    = {}
        options.save_fig    (1,1) logical = true
        options.out_dir     (1,:) char    = pwd
    end

    %% ── Bands to analyse ────────────────────────────────────────────────
    all_fields = fieldnames(struct_CP);
    if isempty(options.bands)
        bands = all_fields;
    else
        bands = options.bands(:);
        for b = 1:numel(bands)
            if ~ismember(bands{b}, all_fields)
                error('Band "%s" not found. Available: %s', ...
                    bands{b}, strjoin(all_fields,', '));
            end
        end
    end
    n_bands = numel(bands);

    fprintf('\n════════════════════════════════════════\n');
    fprintf(' dx=%d | Bands: %s\n', dx, strjoin(bands,', '));
    fprintf('════════════════════════════════════════\n');

    %% ── Create output folder immediately ────────────────────────────────
    if ~isfolder(options.out_dir)
        mkdir(options.out_dir);
        fprintf('Folder created: %s\n', options.out_dir);
    end

    %% ── Colours and styles ───────────────────────────────────────────────
    clr_CP  = [0.78 0.18 0.18];   % red   – infected
    clr_FOL = [0.13 0.55 0.34];   % green – healthy
    fnt     = 'Helvetica';

    %% ── Subplot layout: compact grid with no empty cells ─────────────────
    n_cols = ceil(sqrt(n_bands));
    n_rows = ceil(n_bands / n_cols);
    while n_cols > 1 && ceil(n_bands / (n_cols - 1)) == n_rows
        n_cols = n_cols - 1;
    end
    n_rows = ceil(n_bands / n_cols);

    fprintf('Subplot layout: %d rows x %d columns (for %d bands)\n', ...
            n_rows, n_cols, n_bands);

    fig_CECP = figure('Color','w','Name','CECP – all bands', ...
        'Units','centimeters','Position',[1 1 n_cols*11 n_rows*10]);
    fig_ROC  = figure('Color','w','Name','ROC – all bands', ...
        'Units','centimeters','Position',[1 1 n_cols*11 n_rows*10]);
    fig_STAT = figure('Color','w','Name','Statistics – all bands', ...
        'Units','centimeters','Position',[1 1 n_cols*9  n_rows*10]);

    results = struct();

    %% ── Band loop ────────────────────────────────────────────────────────
    for b = 1:n_bands
        band = bands{b};
        fprintf('\n── Band: %s ──\n', band);

        tbl_CP  = struct_CP.(band);
        tbl_FOL = struct_FOL.(band);
        n_pix_CP  = width(tbl_CP)  - 1;
        n_pix_FOL = width(tbl_FOL) - 1;

        %% ── Compute H and C ──────────────────────────────────────────────
        H_CP  = NaN(n_pix_CP,  1);
        C_CP  = NaN(n_pix_CP,  1);
        H_FOL = NaN(n_pix_FOL, 1);
        C_FOL = NaN(n_pix_FOL, 1);

        fprintf('  CP  (%d pixels)...', n_pix_CP);
        for p = 1:n_pix_CP
            try
                [H_CP(p), C_CP(p)] = complexity_entropy_universal( ...
                    tbl_CP{:,p+1}, dx, options.dy, options.taux, options.tauy, ...
                    options.normalize, options.use_miller);
            catch, end
        end
        fprintf(' done\n');

        fprintf('  FOL (%d pixels)...', n_pix_FOL);
        for p = 1:n_pix_FOL
            try
                [H_FOL(p), C_FOL(p)] = complexity_entropy_universal( ...
                    tbl_FOL{:,p+1}, dx, options.dy, options.taux, options.tauy, ...
                    options.normalize, options.use_miller);
            catch, end
        end
        fprintf(' done\n');

        %% ── Remove NaNs ──────────────────────────────────────────────────
        vCP  = ~isnan(H_CP)  & ~isnan(C_CP);
        vFOL = ~isnan(H_FOL) & ~isnan(C_FOL);
        H_CP  = H_CP(vCP);   C_CP  = C_CP(vCP);
        H_FOL = H_FOL(vFOL); C_FOL = C_FOL(vFOL);

        fprintf('  Valid  CP:%d  FOL:%d\n', numel(H_CP), numel(H_FOL));

        %% ── Statistics (Cohen d only) ────────────────────────────────────
        d_H = cohen_d(H_CP, H_FOL);
        d_C = cohen_d(C_CP, C_FOL);

        fprintf('  Cohen d   d_H=%.3f  d_C=%.3f\n', d_H, d_C);

        %% ── ROC — H and C only ───────────────────────────────────────────
        labels = [ones(numel(H_CP),1); zeros(numel(H_FOL),1)];
        roc_H  = compute_roc(labels, [H_CP; H_FOL], 'H');
        roc_C  = compute_roc(labels, [C_CP; C_FOL], 'C');

        fprintf('  AUC  H=%.3f  C=%.3f\n', roc_H.AUC, roc_C.AUC);

        %% ── CECP subplot ─────────────────────────────────────────────────
        figure(fig_CECP);
        ax1 = subplot(n_rows, n_cols, b);
        hold(ax1,'on'); grid(ax1,'on');
        set(ax1,'FontName',fnt,'FontSize',8,'Box','off');

        if options.show_curves
            plot_limit_curves(ax1, dx);
        end

        plot_ellipse(ax1, H_FOL, C_FOL, clr_FOL, 0.15);
        plot_ellipse(ax1, H_CP,  C_CP,  clr_CP,  0.15);

        scatter(ax1, H_FOL, C_FOL, 20, clr_FOL, 'o', ...
                'filled','MarkerFaceAlpha',0.45,'DisplayName','FOL');
        scatter(ax1, H_CP,  C_CP,  20, clr_CP,  's', ...
                'filled','MarkerFaceAlpha',0.45,'DisplayName','CP');

        plot(ax1, mean(H_FOL), mean(C_FOL), 'o','MarkerSize',8, ...
             'MarkerFaceColor',clr_FOL,'MarkerEdgeColor','k', ...
             'LineWidth',1.0,'HandleVisibility','off');
        plot(ax1, mean(H_CP),  mean(C_CP),  's','MarkerSize',8, ...
             'MarkerFaceColor',clr_CP,'MarkerEdgeColor','k', ...
             'LineWidth',1.0,'HandleVisibility','off');

        xlabel(ax1,'H','FontSize',9);
        ylabel(ax1,'C','FontSize',9);
        title(ax1, band, 'FontSize',9,'FontWeight','bold');
        if b == 1
            legend(ax1,'Location','best','FontSize',7,'Box','off');
        end

        %% ── ROC subplot ──────────────────────────────────────────────────
        figure(fig_ROC);
        ax2 = subplot(n_rows, n_cols, b);
        hold(ax2,'on'); grid(ax2,'on');
        set(ax2,'FontName',fnt,'FontSize',8,'Box','off');

        plot(ax2, roc_H.fpr, roc_H.tpr, 'b-','LineWidth',1.6, ...
             'DisplayName', sprintf('H (AUC=%.3f)', roc_H.AUC));
        plot(ax2, roc_C.fpr, roc_C.tpr, 'r-','LineWidth',1.6, ...
             'DisplayName', sprintf('C (AUC=%.3f)', roc_C.AUC));
        plot(ax2, [0 1],[0 1],'k--','LineWidth',0.7,'HandleVisibility','off');

        xlabel(ax2,'FPR','FontSize',9);
        ylabel(ax2,'TPR','FontSize',9);
        title(ax2, band,'FontSize',9,'FontWeight','bold');
        legend(ax2,'Location','southeast','FontSize',7,'Box','off');

        %% ── Statistics subplot ───────────────────────────────────────────
        figure(fig_STAT);

        % ── Boxplot H ─────────────────────────────────────────────────────
        ax3 = subplot(n_rows, n_cols*2, (b-1)*2+1);
        hold(ax3,'on'); grid(ax3,'on');
        set(ax3,'FontName',fnt,'FontSize',8,'Box','off');

        boxplot(ax3, [H_CP; H_FOL], ...
            [repmat({'CP'},numel(H_CP),1); repmat({'FOL'},numel(H_FOL),1)], ...
            'Colors',[clr_CP; clr_FOL],'Symbol','o','Widths',0.5);
        ylabel(ax3,'H','FontSize',9);
        title(ax3, sprintf('%s – H', band), 'FontSize',8);

        all_H  = [H_CP; H_FOL];
        y_lo_H = prctile(all_H,  2);
        y_hi_H = prctile(all_H, 98);
        pad_H  = max(0.02, 0.08 * (y_hi_H - y_lo_H));
        ylim(ax3, [y_lo_H - pad_H, y_hi_H + pad_H]);

        % ── Boxplot C ─────────────────────────────────────────────────────
        ax4 = subplot(n_rows, n_cols*2, (b-1)*2+2);
        hold(ax4,'on'); grid(ax4,'on');
        set(ax4,'FontName',fnt,'FontSize',8,'Box','off');

        boxplot(ax4, [C_CP; C_FOL], ...
            [repmat({'CP'},numel(C_CP),1); repmat({'FOL'},numel(C_FOL),1)], ...
            'Colors',[clr_CP; clr_FOL],'Symbol','o','Widths',0.5);
        ylabel(ax4,'C','FontSize',9);
        title(ax4, sprintf('%s – C', band), 'FontSize',8);

        all_C  = [C_CP; C_FOL];
        y_lo_C = prctile(all_C,  2);
        y_hi_C = prctile(all_C, 98);
        pad_C  = max(0.02, 0.08 * (y_hi_C - y_lo_C));
        ylim(ax4, [y_lo_C - pad_C, y_hi_C + pad_C]);

        %% ── Save band results ────────────────────────────────────────────
        results.(band).H_CP  = H_CP;
        results.(band).C_CP  = C_CP;
        results.(band).H_FOL = H_FOL;
        results.(band).C_FOL = C_FOL;
        results.(band).stats = struct( ...
            'd_H',  d_H,  'd_C', d_C, ...
            'mean_H_CP',  mean(H_CP),  'std_H_CP',  std(H_CP), ...
            'mean_C_CP',  mean(C_CP),  'std_C_CP',  std(C_CP), ...
            'mean_H_FOL', mean(H_FOL), 'std_H_FOL', std(H_FOL), ...
            'mean_C_FOL', mean(C_FOL), 'std_C_FOL', std(C_FOL));
        results.(band).roc = struct('H', roc_H, 'C', roc_C);
    end

    %% ── Centre last row if incomplete ────────────────────────────────────
    remainder = mod(n_bands, n_cols);
    if remainder > 0
        n_empty   = n_cols - remainder;
        half_shift = n_empty / (2 * n_cols);

        for fig_h = {fig_CECP, fig_ROC}
            fh    = fig_h{1};
            axs   = findobj(fh, 'Type', 'axes');
            for k = 1:remainder
                ax_k = axs(k);
                pos  = get(ax_k, 'Position');
                pos(1) = pos(1) + half_shift;
                set(ax_k, 'Position', pos);
            end
        end

        axs_stat = findobj(fig_STAT, 'Type', 'axes');
        for k = 1 : remainder*2
            ax_k = axs_stat(k);
            pos  = get(ax_k, 'Position');
            pos(1) = pos(1) + half_shift;
            set(ax_k, 'Position', pos);
        end
    end

    %% ── Summary table ────────────────────────────────────────────────────
    T = build_summary_table(results, bands);

    % Print to terminal — comparison columns only
    fprintf('\n════════════════════════════════════════════════════════════════\n');
    fprintf(' CECP Summary  –  dx = %d\n', dx);
    fprintf('════════════════════════════════════════════════════════════════\n');
    fprintf('%-8s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s\n', ...
            'Band','d_H','AUC_H','thr_H','TPR_H','FPR_H', ...
                   'd_C','AUC_C','thr_C','TPR_C','FPR_C');
    fprintf('%-8s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s\n', ...
            '--------','------','------','------','------','------', ...
                       '------','------','------','------','------');
    for r = 1:height(T)
        fprintf('%-8s  %6.3f  %6.3f  %6.4f  %6.3f  %6.3f  %6.3f  %6.3f  %6.4f  %6.3f  %6.3f\n', ...
            T.Band{r}, ...
            T.d_H(r),   T.AUC_H(r), T.thr_H(r), T.TPR_H(r), T.FPR_H(r), ...
            T.d_C(r),   T.AUC_C(r), T.thr_C(r), T.TPR_C(r), T.FPR_C(r));
    end
    fprintf('════════════════════════════════════════════════════════════════\n\n');

    % Graphical display
    [fig_TAB_H, fig_TAB_C] = show_summary_table(T, dx);

    %% ── Save figures and table ───────────────────────────────────────────
    if options.save_fig
        base = fullfile(options.out_dir, sprintf('CECP_allbands_dx%d', dx));

        figs  = {fig_CECP, fig_ROC,  fig_STAT,  fig_TAB_H,   fig_TAB_C};
        names = {'_CECP',  '_ROC',   '_STAT',   '_TABLE_H',  '_TABLE_C'};

        for fi = 1:numel(figs)
            exportgraphics(figs{fi}, [base names{fi} '.png'], ...
                           'Resolution',300,'BackgroundColor','white');
            exportgraphics(figs{fi}, [base names{fi} '.pdf'], ...
                           'ContentType','vector','BackgroundColor','white');
            savefig(figs{fi}, [base names{fi} '.fig']);
        end

        save([base '_results.mat'], 'results');

        fprintf('\nFigures saved in: %s\n', options.out_dir);
        for fi = 1:numel(figs)
            fprintf('  %s%s.{png,pdf,fig}\n', base, names{fi});
        end
        fprintf('  %s_results.mat\n', base);
        fprintf('  %s_summary.csv\n', base);
    end

    % Save table to CSV and into results struct
    if options.save_fig
        csv_path = fullfile(options.out_dir, sprintf('CECP_allbands_dx%d_summary.csv', dx));
        writetable(T, csv_path);
    end
    results.summary = T;
end


%% ════════════════════════════════════════════════════════════════════════
%%  HELPERS
%% ════════════════════════════════════════════════════════════════════════

function r = compute_roc(labels, scores, name)
    [fpr, tpr, thr, auc] = perfcurve(labels, scores, 1);
    if auc < 0.5
        [fpr, tpr, thr, auc] = perfcurve(labels, -scores, 1);
        r.inverted  = true;
        r.direction = sprintf('%s: CP < FOL', name);
    else
        r.inverted  = false;
        r.direction = sprintf('%s: CP > FOL', name);
    end
    % Optimal threshold: Youden criterion (max TPR - FPR)
    % perfcurve adds a sentinel row with thr=Inf → exclude it
    valid      = isfinite(thr);
    [~, idx]   = max(tpr(valid) - fpr(valid));
    idx_valid  = find(valid);
    opt_idx    = idx_valid(idx);
    r.fpr      = fpr;
    r.tpr      = tpr;
    r.thresh   = thr;
    r.AUC      = auc;
    r.opt_thr  = thr(opt_idx);
    r.opt_tpr  = tpr(opt_idx);
    r.opt_fpr  = fpr(opt_idx);
end

function d = cohen_d(x, y)
    nx = numel(x); ny = numel(y);
    sp = sqrt(((nx-1)*var(x) + (ny-1)*var(y)) / (nx+ny-2));
    d  = (mean(x) - mean(y)) / sp;
end

function plot_ellipse(ax, H, C, col, alpha_val)
    if numel(H) < 3, return; end
    mu = [mean(H), mean(C)];
    S  = cov(H, C);
    [V, D] = eig(S);
    theta  = linspace(0, 2*pi, 300);
    ell    = V * (2*sqrt(D) * [cos(theta); sin(theta)]);
    fill(ax, mu(1)+ell(1,:), mu(2)+ell(2,:), col, ...
         'FaceAlpha',alpha_val,'EdgeColor',col, ...
         'LineWidth',1.2,'HandleVisibility','off');
end

function s = fmt_p(p)
% FMT_P  Format a p-value: decimal if p>=0.001, scientific notation otherwise.
    if p >= 0.001
        s = sprintf('%.3f', p);
    else
        s = sprintf('%.2e', p);
        s = regexprep(s, 'e([+-])0*(\d+)$', 'e$1$2');
    end
end

function plot_limit_curves(ax, dx)
    try
        c_max = maximum_complexity_entropy(dx, 1, 200, false);
        c_min = minimum_complexity_entropy(dx, 1, 200, false);
        plot(ax, c_max(:,1), c_max(:,2),'-','Color',[0.6 0.6 0.6], ...
             'LineWidth',1.2,'HandleVisibility','off');
        plot(ax, c_min(:,1), c_min(:,2),'--','Color',[0.6 0.6 0.6], ...
             'LineWidth',1.2,'HandleVisibility','off');
    catch
        warning('Limit curves not available.');
    end
end

function T = build_summary_table(results, bands)
% BUILD_SUMMARY_TABLE  Assembles a summary table for all bands.
%
% Columns:
%   Band           : band name
%   mean_H_CP/FOL  : mean H for CP and FOL
%   std_H_CP/FOL   : std  H for CP and FOL
%   d_H            : Cohen d for H
%   AUC_H          : ROC AUC for H
%   thr_H          : optimal threshold H (Youden)
%   TPR_H          : TPR at optimal threshold H
%   FPR_H          : FPR at optimal threshold H
%   (same columns for C)

    n = numel(bands);

    Band       = bands(:);
    mean_H_CP  = zeros(n,1);  std_H_CP  = zeros(n,1);
    mean_H_FOL = zeros(n,1);  std_H_FOL = zeros(n,1);
    d_H        = zeros(n,1);  AUC_H     = zeros(n,1);
    thr_H      = zeros(n,1);  TPR_H     = zeros(n,1);  FPR_H = zeros(n,1);
    mean_C_CP  = zeros(n,1);  std_C_CP  = zeros(n,1);
    mean_C_FOL = zeros(n,1);  std_C_FOL = zeros(n,1);
    d_C        = zeros(n,1);  AUC_C     = zeros(n,1);
    thr_C      = zeros(n,1);  TPR_C     = zeros(n,1);  FPR_C = zeros(n,1);

    for b = 1:n
        bnd = bands{b};
        s   = results.(bnd).stats;
        rH  = results.(bnd).roc.H;
        rC  = results.(bnd).roc.C;

        mean_H_CP(b)  = s.mean_H_CP;   std_H_CP(b)  = s.std_H_CP;
        mean_H_FOL(b) = s.mean_H_FOL;  std_H_FOL(b) = s.std_H_FOL;
        d_H(b)        = s.d_H;
        AUC_H(b)      = rH.AUC;
        thr_H(b)      = rH.opt_thr;
        TPR_H(b)      = rH.opt_tpr;
        FPR_H(b)      = rH.opt_fpr;

        mean_C_CP(b)  = s.mean_C_CP;   std_C_CP(b)  = s.std_C_CP;
        mean_C_FOL(b) = s.mean_C_FOL;  std_C_FOL(b) = s.std_C_FOL;
        d_C(b)        = s.d_C;
        AUC_C(b)      = rC.AUC;
        thr_C(b)      = rC.opt_thr;
        TPR_C(b)      = rC.opt_tpr;
        FPR_C(b)      = rC.opt_fpr;
    end

    T = table(Band, ...
        mean_H_CP, std_H_CP, mean_H_FOL, std_H_FOL, d_H, AUC_H, thr_H, TPR_H, FPR_H, ...
        mean_C_CP, std_C_CP, mean_C_FOL, std_C_FOL, d_C, AUC_C, thr_C, TPR_C, FPR_C, ...
        'VariableNames', { ...
            'Band', ...
            'mean_H_CP','std_H_CP','mean_H_FOL','std_H_FOL','d_H','AUC_H','thr_H','TPR_H','FPR_H', ...
            'mean_C_CP','std_C_CP','mean_C_FOL','std_C_FOL','d_C','AUC_C','thr_C','TPR_C','FPR_C'});
end

function [fig_H, fig_C] = show_summary_table(T, dx)
% SHOW_SUMMARY_TABLE  Creates two separate figures, one for H and one for C.
%
%   fig_H:  Band | d_H | AUC_H | thr_H | TPR_H | FPR_H
%   fig_C:  Band | d_C | AUC_C | thr_C | TPR_C | FPR_C

    n   = height(T);
    fnt = 'Helvetica';

    cols_H = { ...
        'Band',  @(r) T.Band{r};                    ...
        'd_H',   @(r) sprintf('%.3f', T.d_H(r));   ...
        'AUC_H', @(r) sprintf('%.3f', T.AUC_H(r)); ...
        'thr_H', @(r) sprintf('%.4f', T.thr_H(r)); ...
        'TPR_H', @(r) sprintf('%.3f', T.TPR_H(r)); ...
        'FPR_H', @(r) sprintf('%.3f', T.FPR_H(r)); ...
    };
    cols_C = { ...
        'Band',  @(r) T.Band{r};                    ...
        'd_C',   @(r) sprintf('%.3f', T.d_C(r));   ...
        'AUC_C', @(r) sprintf('%.3f', T.AUC_C(r)); ...
        'thr_C', @(r) sprintf('%.4f', T.thr_C(r)); ...
        'TPR_C', @(r) sprintf('%.3f', T.TPR_C(r)); ...
        'FPR_C', @(r) sprintf('%.3f', T.FPR_C(r)); ...
    };

    fig_H = make_table_fig(cols_H, n, fnt, sprintf('H  –  dx = %d', dx));
    fig_C = make_table_fig(cols_C, n, fnt, sprintf('C  –  dx = %d', dx));
end

function fig = make_table_fig(cols, n, fnt, title_str)
% MAKE_TABLE_FIG  Creates a single figure with the specified table.

    nc = size(cols, 1);

    % Column widths: Band wider, metrics distributed evenly
    w_band   = 0.16;
    w_metric = (0.96 - w_band - 0.008*(nc-1)) / (nc-1);
    col_w    = [w_band, repmat(w_metric, 1, nc-1)];
    gap_col  = 0.008;
    x_margin = 0.02;

    % Adaptive cell height
    cell_h = min(0.09, 0.80 / (n + 1));

    % Figure dimensions
    fig_w = 16;
    fig_h = max(6, (n + 2) * 0.65 + 1.5);

    fig = figure('Color','w', ...
        'Name', title_str, ...
        'Units','centimeters', ...
        'Position',[2 2 fig_w fig_h]);
    ax = axes(fig,'Position',[0 0 1 1],'Visible','off');
    ax.XLim = [0 1]; ax.YLim = [0 1];
    hold(ax,'on');

    % Compute x_starts
    x_starts = zeros(1, nc);
    xl = x_margin;
    for c = 1:nc
        x_starts(c) = xl;
        xl = xl + col_w(c) + gap_col;
    end

    y_top = 0.88;

    % Header
    clr_head = [0.88 0.88 0.88];
    clr_line = [0.40 0.40 0.40];
    y_h = y_top - cell_h;
    for c = 1:nc
        w = col_w(c);
        rectangle(ax,'Position',[x_starts(c), y_h, w, cell_h], ...
                  'FaceColor',clr_head,'EdgeColor',clr_line,'LineWidth',0.7);
        text(ax, x_starts(c)+w/2, y_h+cell_h/2, cols{c,1}, ...
             'HorizontalAlignment','center','VerticalAlignment','middle', ...
             'FontName',fnt,'FontSize',8,'FontWeight','bold','Color',[0 0 0]);
    end

    % Data rows
    for r = 1:n
        y0 = y_top - (r+1)*cell_h;
        for c = 1:nc
            w = col_w(c);
            rectangle(ax,'Position',[x_starts(c), y0, w, cell_h], ...
                      'FaceColor',[1 1 1],'EdgeColor',clr_line,'LineWidth',0.5);
            fw  = pick_ternary(c == 1, 'bold', 'normal');
            txt = cols{c,2}(r);
            text(ax, x_starts(c)+w/2, y0+cell_h/2, txt, ...
                 'HorizontalAlignment','center','VerticalAlignment','middle', ...
                 'FontName',fnt,'FontSize',8,'FontWeight',fw,'Color',[0 0 0]);
        end
    end

    title(ax, title_str, ...
          'FontName',fnt,'FontSize',10,'FontWeight','bold', ...
          'Units','normalized','Position',[0.5 0.96 0]);
end

% ── Ternary operator ─────────────────────────────────────────────────────
function v = pick_ternary(cond, a, b)
    if cond, v = a; else, v = b; end
end
