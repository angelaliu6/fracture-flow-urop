outDir = 'transport';
if ~exist(outDir,'dir'); mkdir(outDir); end

Lx = 150; Ly = 50;
visc = 1e-3; pinlet = 1e5; poutlet = 0; tol = 1e-12; seed = NaN;

num_param_sets = 100;
num_trials     = 1000;  % trials per param set (each trial runs Np particles)
Np             = 1000;  % particles per trial
timeout        = 1e100;

weights = rand(num_param_sets, 3);

all_trial_bkt_data = cell(num_param_sets, 1);  % {paramset}{trial} = sorted finite BKTs

outpdf_cdf = fullfile(outDir, 'bkt_1-CDF_all_all.pdf');
if exist(outpdf_cdf, 'file'); delete(outpdf_cdf); end

for paramset_idx = 1:num_param_sets

    discspacing  = 1   + (3   - 1  ) * weights(paramset_idx, 1); % [1, 3]
    length_ratio = 0.1 + (0.5 - 0.1) * weights(paramset_idx, 2); % [0.1, 0.5]
    b            = 0.01e-3 + (0.3e-3 - 0.01e-3) * weights(paramset_idx, 3); % [0.01e-3, 0.3e-3]

    lx = length_ratio * Lx;
    ly = length_ratio * Ly;
    Nx = round(Lx / discspacing);
    Ny = round(Ly / discspacing);
    bmin = b;
    bmax = b;

    trial_bkt_data = cell(num_trials, 1);
    connected_count = 0;

    for trial_idx = 1:num_trials

        [linktrans,linkapert,linknods,backbonenods,isolatednods, ...
                conninletnods,connoutletnods,pvec,tol,linkflow,linkvel, ...
                conninletlinks,flow_in] = ...
            run_trial_detailed(Lx,Ly,Nx,Ny,lx,ly,bmin,bmax,visc,pinlet,poutlet,tol,seed);

        if flow_in > 0
            connected_count = connected_count + 1;
            % posnode is injected into this workspace by run_trial_detailed via assignin
            [npaths, tpaths] = simul_transport(linkflow,linkvel,linknods,posnode, ...
                conninletnods,conninletlinks,connoutletnods,tol,Np,timeout);
            breakthrough_times = cellfun(@(t) t(end), tpaths);
            trial_bkt_data{trial_idx} = sort(breakthrough_times(isfinite(breakthrough_times)));
        end

    end

    p0 = 1 - connected_count / num_trials;
    all_trial_bkt_data{paramset_idx} = struct( ...
        'trial_data', {trial_bkt_data}, ...
        'p0', p0);

    % one plot per param set: num_trials curves, each curve = 1-CDF of Np particles
    fig_cdf = figure(1); clf;
    hold on;
    for trial_idx = 1:num_trials
        fb = trial_bkt_data{trial_idx};
        if ~isempty(fb)
            n_fb = numel(fb);
            survival = (n_fb - (1:n_fb)') / n_fb;
            loglog(fb, survival, 'Color', [0.3 0.3 0.8 0.3]);
        end
    end
    hold off;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    ax = gca;
    ax.XAxis.TickLabelFormat = '10^{%g}';
    xlabel('Breakthrough time (s)');
    ylabel('1 - CDF');
    title(sprintf('Param set %d  (disc=%.2f, lr=%.2f, b=%.2e, p0=%.2f)', ...
        paramset_idx, discspacing, length_ratio, b, p0));
    grid on;
    exportgraphics(fig_cdf, outpdf_cdf, 'Append', true);

end

save(fullfile(outDir, 'all_trial_bkt_data.mat'), 'all_trial_bkt_data', 'weights', '-v7.3');
