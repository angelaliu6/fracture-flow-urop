outDir = '3-parameter data';
if ~exist(outDir,'dir'); mkdir(outDir); end

Lx = 150;
Ly = 50;
visc = 1e-3;
pinlet = 1e5;
poutlet = 0;
tol = 1e-12;
seed = NaN;
num_trials = 1000; % for each (discspacing, lx, ly, b) --> 1000
num_param_sets = 10000 % --> 100
weights = rand(num_param_sets, 3) % 3 = number of parameters -- can be changed to only be in the for loop

parfor paramset_idx = 1:num_param_sets % IF YOU ARE ADDING MORE PARAMETER SETS, SHIFT THIS INDEX
    discspacing = 5*weights(paramset_idx,1)+0.1*(1-weights(paramset_idx,1))
    b = 0.01e-3*weights(paramset_idx,2)+0.3e-3*(1-weights(paramset_idx,2))
    % lx = 60*weights(paramset_idx,2)+80*(1-weights(paramset_idx,3))
    % ly = 20*weights(paramset_idx,3)+30*(1-weights(paramset_idx,4))
    length_ratio = 0.02*weights(paramset_idx, 3)+0.5*(1-weights(paramset_idx, 3))
    lx = length_ratio*Lx
    ly = length_ratio*Ly
    Nx = round(Lx/discspacing)
    Ny = round(Ly/discspacing)
   
    bmin = b;
    bmax = b;

    % logflows = zeros(num_trials, 1)
    flows = zeros(num_trials,1)
    for trial_idx = 1:num_trials
        [is_connected, flow] = run_trial(Lx, Ly, Nx, Ny, lx, ly, bmin, bmax, visc, pinlet, poutlet, tol, seed)
        % if is_connected
        %     logflows(trial_idx) = log(flow)
        % else
        %     logflows(trial_idx) = -999
        % end
        flows(trial_idx) = flow
    end
    
    tables{paramset_idx} = table( ...
        repmat(discspacing, num_trials, 1), ...
        repmat(length_ratio,num_trials, 1), ...
        repmat(b,           num_trials, 1), ...
        flows, ...
        'VariableNames', {'discspacing','length_ratio','b','Q'} ...
    );
    filename = "flow_outputs_" + paramset_idx + ".csv";
    writetable(tables{paramset_idx}, filename);
end
