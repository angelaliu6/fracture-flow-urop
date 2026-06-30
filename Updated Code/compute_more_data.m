Lx = 150;
Ly = 50;
visc = 1e-3;
pinlet = 1e5;
poutlet = 0;
tol = 1e-12;
seed = NaN;
num_trials = 1000;
% redo the files for 'bad' parameter sets with p0>0.9, using 10,000 trials instead
get_more_data_for = readmatrix(fullfile(pwd, '3-parameter data', 'get_more_data_for.csv'));
num_param_sets = size(get_more_data_for, 1);
parfor paramset_idx = 1:num_param_sets
    file_idx = get_more_data_for(paramset_idx, 1)
    discspacing = get_more_data_for(paramset_idx, 2)
    length_ratio = get_more_data_for(paramset_idx, 3)
    b = get_more_data_for(paramset_idx, 4)
    bmin = b;
    bmax = b;
    lx = length_ratio*Lx
    ly = length_ratio*Ly
    Nx = round(Lx/discspacing)
    Ny = round(Ly/discspacing)

    p0 = get_more_data_for(paramset_idx, 5)
    larger_num_trials = floor(num_trials/5 / (1-p0+0.01))
    flows = zeros(larger_num_trials,1)
    for trial_idx = 1:larger_num_trials
        [is_connected, flow] = run_trial(Lx, Ly, Nx, Ny, lx, ly, bmin, bmax, visc, pinlet, poutlet, tol, seed)
        flows(trial_idx) = flow
    end
    
    tables{paramset_idx} = table( ...
        repmat(discspacing, larger_num_trials, 1), ...
        repmat(length_ratio,larger_num_trials, 1), ...
        repmat(b,           larger_num_trials, 1), ...
        flows, ...
        'VariableNames', {'discspacing','length_ratio','b','Q'} ...
    );
filename = fullfile(pwd, '3-parameter data', sprintf('flow_outputs_%d.csv', file_idx));
end
