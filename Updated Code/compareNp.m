% Compare 1000 vs 10000 particles on the same network.
% Uses first 5 param sets from existing weights; fixed seeds guarantee the
% same network realization is used for both particle counts.
% Outputs: transport/compare_Np.mat

outDir = 'transport';
tmp    = load(fullfile(outDir, 'all_trial_bkt_data.mat'), 'weights');
weights = tmp.weights;

Lx = 150; Ly = 50;
visc = 1e-3; pinlet = 1e5; poutlet = 0; tol = 1e-12; timeout = 1e100;

N_PS    = 5;
N_TRIAL = 5;

bkt_1000  = cell(N_PS, N_TRIAL);
bkt_10000 = cell(N_PS, N_TRIAL);
connected = false(N_PS, N_TRIAL);

for ps_idx = 1:N_PS

    discspacing  = 1.0     + (5.0     - 1.0    ) * weights(ps_idx, 1);
    length_ratio = 0.3     + (0.8     - 0.3    ) * weights(ps_idx, 2);
    b            = 0.01e-3 + (0.3e-3  - 0.01e-3) * weights(ps_idx, 3);

    lx = length_ratio * Lx;
    ly = length_ratio * Ly;
    Nx = round(Lx / discspacing);
    Ny = round(Ly / discspacing);
    bmin = b; bmax = b;

    for tr_idx = 1:N_TRIAL

        seed = 100 * ps_idx + tr_idx;   % fixed seed → same network for both Np

        [linktrans, linkapert, linknods, backbonenods, isolatednods, ...
                conninletnods, connoutletnods, pvec, tol, linkflow, linkvel, ...
                conninletlinks, flow_in] = ...
            run_trial_detailed(Lx, Ly, Nx, Ny, lx, ly, bmin, bmax, ...
                               visc, pinlet, poutlet, tol, seed);

        if flow_in > 0
            connected(ps_idx, tr_idx) = true;

            for Np = [1000, 10000]
                [~, tpaths] = simul_transport(linkflow, linkvel, linknods, posnode, ...
                    conninletnods, conninletlinks, connoutletnods, tol, Np, timeout);
                bt = sort(cellfun(@(t) t(end), tpaths));
                bt = bt(isfinite(bt));
                if Np == 1000
                    bkt_1000{ps_idx, tr_idx}  = bt;
                else
                    bkt_10000{ps_idx, tr_idx} = bt;
                end
            end

            fprintf('PS %d, trial %d: connected  (Np=1000: %d,  Np=10000: %d)\n', ...
                ps_idx, tr_idx, ...
                numel(bkt_1000{ps_idx, tr_idx}), numel(bkt_10000{ps_idx, tr_idx}));
        else
            fprintf('PS %d, trial %d: disconnected\n', ps_idx, tr_idx);
        end

    end
end

save(fullfile(outDir, 'compare_Np.mat'), ...
     'bkt_1000', 'bkt_10000', 'connected', 'weights', '-v7.3');
fprintf('Saved compare_Np.mat\n');
