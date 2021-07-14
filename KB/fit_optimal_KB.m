function fit_optimal_KB(results_dir)

    load(sprintf('%s/KB/Hateren.mat', results_dir), 'S7_unique');
    [C, alpha] = get_DCT_proj_matrix(true, 10, 10, 20, 10); 
    N_trials = 100;
    
    for i = 1:N_trials        
        if (i == 1)
            [theta_est{i}, phi_est{i}, v_norm{i}] = get_Klein_coords_opt(S7_unique, C, alpha(i,:));
        else
            [theta_est{i}, phi_est{i}, v_norm{i}] = get_Klein_coords_opt(S7_unique, C, alpha(i,:), theta_est{i-1}, phi_est{i-1});
        end
        [alpha(i+1,:), ~, v_err{i}, ~, v_norm{i+1}] = fit_Klein_bottle(S7_unique, C, theta_est{i}, phi_est{i}, v_norm{i});
        if (i>1 && (sum(v_err{i})>sum(v_err{i-1})))
            break
        end
    end

    save(sprintf('%s/KB/Fit_1_10_10_20_10.mat', results_dir), 'C', 'alpha', 'N_trials', 'theta_est', 'phi_est', 'v_norm', 'v_err', '-v7.3');

end