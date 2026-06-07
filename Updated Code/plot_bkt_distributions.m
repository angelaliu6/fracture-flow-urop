outDir = 'transport';
num_param_sets = 100;
stats  = {'min_bkt', 'median_bkt', 'max_bkt'};
colors = {[0.2 0.4 0.8], [0.9 0.4 0.1], [0.2 0.7 0.3]};

% first pass: read all CSVs to find global log10 x-axis limits
all_finite = [];
for n = 1:num_param_sets
    fname = fullfile(outDir, sprintf('transport_outputs_%d.csv', n));
    if ~exist(fname, 'file'); continue; end
    T = readtable(fname);
    for s = 1:length(stats)
        vals = T.(stats{s});
        all_finite = [all_finite; vals(isfinite(vals) & vals > 0)];
    end
end

if isempty(all_finite)
    error('No finite breakthrough times found in any CSV.');
end
xmin_log = floor(log10(min(all_finite)));
xmax_log = ceil(log10(max(all_finite)));
tick_vals = xmin_log:xmax_log;

% second pass: one overlaid plot per param set, appended into a single PDF
outpdf = fullfile(outDir, 'bkt_distributions.pdf');
if exist(outpdf, 'file'); delete(outpdf); end  % so Append starts fresh

for n = 1:num_param_sets
    fname = fullfile(outDir, sprintf('transport_outputs_%d.csv', n));
    if ~exist(fname, 'file'); continue; end
    T = readtable(fname);

    fig = figure(1); clf;
    hold on;
    for s = 1:length(stats)
        col = stats{s};
        vals = T.(col);
        finite_v = vals(isfinite(vals) & vals > 0);
        if ~isempty(finite_v)
            histogram(log10(finite_v), ...
                'BinLimits', [xmin_log, xmax_log], 'BinWidth', 0.1, ...
                'FaceColor', colors{s}, 'FaceAlpha', 0.5, ...
                'DisplayName', strrep(col, '_', ' '));
        end
    end
    hold off;
    xlim([xmin_log, xmax_log]);
    set(gca, 'XTick', tick_vals, ...
        'XTickLabel', arrayfun(@(x) sprintf('10^{%d}', x), tick_vals, 'UniformOutput', false));
    xlabel('Breakthrough time (s)');
    ylabel('Count');
    legend('Location', 'best');
    title(sprintf('Param set %d  (disc=%.2f, lr=%.2f, b=%.2e)', ...
        n, T.discspacing(1), T.length_ratio(1), T.b(1)));

    exportgraphics(fig, outpdf, 'Append', true);
end

fprintf('Saved %s\n', outpdf);
