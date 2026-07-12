outDir = 'transport modified param ranges';
if ~exist(outDir,'dir'); mkdir(outDir); end

Lx = 150; Ly = 50;
visc = 1e-3; pinlet = 1e5; poutlet = 0; tol = 1e-12; seed = NaN;

num_param_sets  = 5;
num_trials_max  = 10000;  % run once per param set; first 1000 trials double as the num_trials=1000 sample
Np              = 1000;   % particles per trial
timeout         = 1e100;

weights = rand(num_param_sets, 3);

all_trial_bkt_data_10000 = cell(num_param_sets, 1);  % {paramset}{trial} = sorted finite BKTs, 10000 trials

for paramset_idx = 1:num_param_sets

    discspacing  = 1*weights(paramset_idx, 1) + 3*(1 - weights(paramset_idx, 1));    % [1, 3]
    b            = 0.01e-3*weights(paramset_idx, 2) + 0.3e-3*(1 - weights(paramset_idx, 2)); % [0.01e-3, 0.3e-3]
    length_ratio = 0.02*weights(paramset_idx, 3) + 0.5*(1 - weights(paramset_idx, 3));  % [0.02, 0.5]

    lx = length_ratio * Lx;
    ly = length_ratio * Ly;
    Nx = round(Lx / discspacing);
    Ny = round(Ly / discspacing);
    bmin = b;
    bmax = b;

    trial_bkt_data = cell(num_trials_max, 1);
    connected_count = 0;

    for trial_idx = 1:num_trials_max

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

        if mod(trial_idx, 1000) == 0
            fprintf('param set %d: %d / %d trials done\n', paramset_idx, trial_idx, num_trials_max);
        end

    end

    p0 = 1 - connected_count / num_trials_max;
    all_trial_bkt_data_10000{paramset_idx} = struct( ...
        'trial_data', {trial_bkt_data}, ...
        'p0', p0);

end

% the num_trials=1000 sample for each param set is just the first 1000 trials
% of the num_trials=10000 run above (same underlying weights/params, nested
% subsample) -- this is reconstructed in the analysis notebook rather than
% duplicated here.
save(fullfile(outDir, 'compare_num_trials.mat'), 'all_trial_bkt_data_10000', 'weights', '-v7.3');
