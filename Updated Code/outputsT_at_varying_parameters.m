outDir = 'transport';
if ~exist(outDir,'dir'); mkdir(outDir); end

Lx = 150; Ly = 50;
visc = 1e-3; pinlet = 1e5; poutlet = 0; tol = 1e-12; seed = NaN;

num_param_sets = 100;
num_trials     = 100;  % trials per param set (each trial runs Np particles)
Np             = 100;  % particles per trial
timeout        = 999999;

weights = rand(num_param_sets, 3);

all_min_bkt    = cell(num_param_sets, 1);
all_median_bkt = cell(num_param_sets, 1);
all_max_bkt    = cell(num_param_sets, 1);

for paramset_idx = 1:num_param_sets

    discspacing  = 1   + (5   - 1  ) * weights(paramset_idx, 1);
    length_ratio = 0.3 + (0.8 - 0.3) * weights(paramset_idx, 2);
    b            = 0.01e-3 + (0.3e-3 - 0.01e-3) * weights(paramset_idx, 3);

    lx = length_ratio * Lx;
    ly = length_ratio * Ly;
    Nx = round(Lx / discspacing);
    Ny = round(Ly / discspacing);
    bmin = b;
    bmax = b;

    min_bkt    = zeros(num_trials, 1);
    median_bkt = zeros(num_trials, 1);
    max_bkt    = zeros(num_trials, 1);

    for trial_idx = 1:num_trials

        [linktrans,linkapert,linknods,backbonenods,isolatednods, ...
                conninletnods,connoutletnods,pvec,tol,linkflow,linkvel, ...
                conninletlinks,flow_in] = ...
            run_trial_detailed(Lx,Ly,Nx,Ny,lx,ly,bmin,bmax,visc,pinlet,poutlet,tol,seed);

        % posnode is injected into this workspace by run_trial_detailed via assignin
        [npaths, tpaths] = simul_transport(linkflow,linkvel,linknods,posnode, ...
            conninletnods,conninletlinks,connoutletnods,tol,Np,timeout);

        breakthrough_times = cellfun(@(t) t(end), tpaths);
        min_bkt(trial_idx)    = min(breakthrough_times);
        median_bkt(trial_idx) = median(breakthrough_times);
        max_bkt(trial_idx)    = max(breakthrough_times);

        % histogram of particle breakthrough times for this specific trial
        finite_bkt = breakthrough_times(isfinite(breakthrough_times));
        if ~isempty(finite_bkt)
            figure(1); clf;
            histogram(log10(finite_bkt), 'NumBins', 20);
            ax = gca;
            ax.XTick = unique(round(ax.XTick));
            ax.XTickLabel = arrayfun(@(x) sprintf('10^{%d}',x), ax.XTick, 'UniformOutput', false);
            xlabel('Breakthrough time (s)');
            title("Param set " + paramset_idx + ", trial " + trial_idx);
            save_figure(1, outDir, "bkt_histogram_" + paramset_idx + "_" + trial_idx);
        end

    end

    all_min_bkt{paramset_idx}    = min_bkt;
    all_median_bkt{paramset_idx} = median_bkt;
    all_max_bkt{paramset_idx}    = max_bkt;

    T = table( ...
        repmat(discspacing,  num_trials, 1), ...
        repmat(length_ratio, num_trials, 1), ...
        repmat(b,            num_trials, 1), ...
        min_bkt, median_bkt, max_bkt, ...
        'VariableNames', {'discspacing','length_ratio','b','min_bkt','median_bkt','max_bkt'} ...
    );
    writetable(T, fullfile(outDir, "transport_outputs_" + paramset_idx + ".csv"));

end

% compute global log10 x-axis limits across all param sets and all three stats
all_vals = [cell2mat(all_min_bkt); cell2mat(all_median_bkt); cell2mat(all_max_bkt)];
finite_vals = all_vals(isfinite(all_vals) & all_vals > 0);
xmin_log = floor(log10(min(finite_vals)));
xmax_log = ceil(log10(max(finite_vals)));
tick_vals = xmin_log:xmax_log;

for paramset_idx = 1:num_param_sets
    fig = figure(2); clf;
    subtitles = {"Min BKT", "Median BKT", "Max BKT"};
    data = {all_min_bkt{paramset_idx}, all_median_bkt{paramset_idx}, all_max_bkt{paramset_idx}};
    for k = 1:3
        subplot(1, 3, k);
        vals = data{k};
        finite_v = vals(isfinite(vals) & vals > 0);
        if ~isempty(finite_v)
            histogram(log10(finite_v), 'BinLimits', [xmin_log, xmax_log], 'NumBins', 20);
        end
        xlim([xmin_log, xmax_log]);
        set(gca, 'XTick', tick_vals, ...
            'XTickLabel', arrayfun(@(x) sprintf('10^{%d}', x), tick_vals, 'UniformOutput', false));
        xlabel('Breakthrough time (s)');
        title(subtitles{k});
    end
    sgtitle("Param set " + paramset_idx);
    set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0.5 0.5 10 3]);
    print(fig, '-dpdf', '-loose', fullfile(outDir, char("bkt_dist_" + paramset_idx + ".pdf")));
end
