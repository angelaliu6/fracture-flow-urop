function [info_file, probability_connected, log_total_flow_rate, mean_log_total_flow_rate, std_log_total_flow_rate] = info_output_at_parameters(paramset_idx, Lx, Ly, Nx, Ny, lx, ly, bmin, bmax, visc, pinlet, poutlet, tol, seed, num_trials, outDir)
    connectedcounter = 0;
    total_flow_rate = zeros(num_trials,1); % will only end up storing the total_flow_rate for connectedcounter many trials
                
    is_connected_trial = false;
    for trial = 1:num_trials
        [is_connected_trial, total_flow_rate_trial] = run_trial(Lx, Ly, Nx, Ny, lx, ly, bmin, bmax, visc, pinlet, poutlet, tol, seed)
        if is_connected_trial
            connectedcounter = connectedcounter + 1;
            total_flow_rate(connectedcounter) = total_flow_rate_trial;
        end
    end
    % total_flow_rate = total_flow_rate(1:connectedcounter)

    probability_connected = connectedcounter / num_trials;
    log_total_flow_rate = log( total_flow_rate(1:connectedcounter) ); %
    mean_log_total_flow_rate = mean(data);
    std_log_total_flow_rate = std(data);
    baseFileName = sprintf('run_data_%d.mat', paramset_idx);
    if probability_connected > 0
        save(fullfile(outDir, baseFileName),'Ny', 'lx', 'ly', 'bmin', 'probability_connected', 'log_total_flow_rate', 'mean_log_total_flow_rate', 'std_log_total_flow_rate')
    end
    [info_file] = baseFileName;

end