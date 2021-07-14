function generate_image_results(results_dir)
   
    se_targ = 0.1;
    
    %% Original image manifold
    
    load(sprintf('%s/KB/Hateren.mat', results_dir), 'S7_unique');    
    rng(404293);
    calib_points = randperm(size(S7_unique,1), 10000)';
    manifold_curvature(sprintf('%s/KB/Hateren_KNN_100_10', results_dir), S7_unique, 2, se_targ, 'calib_points', calib_points);
   
    %% Carlsson k0 embedding and with noise
    
    load(sprintf('%s/KB/HaterenCarlssonMap.mat', results_dir), 'theta_map', 'phi_map');

    load(sprintf('%s/KB/CarlssonKB.mat', results_dir), 'Klein_R9_fcn');
    X = arrayfun(@(t,p) Klein_R9_fcn(t,p), theta_map, phi_map, 'un', false);
    X = vertcat(X{:});
    
    DCT = load(sprintf('%s/KB/DCTParams.mat', results_dir));
    
    median_Enorm = median(sqrt(sum(X.^2,2)));
    
    rng(240240);
    noise = randn(size(X));

    sig = [0.000 0.007 0.010 0.030];
    
    for i = 1:length(sig)
        v = X+sig(i)*median_Enorm*noise;
        v = v*DCT.A*DCT.lambda;
        v = v./sqrt(sum(v.^2,2));
        if (i == 1)
            v = v(:,1:5);
            load(sprintf('%s/KB/Hateren_KNN_100_10/Curvature.mat', results_dir), 'ball_r');
            curvature_at_length_scale(sprintf('%s/KB/Carlsson_KNN_100_10_Data_Radius', results_dir), v, 2, ball_r, calib_points);
        end
        manifold_curvature(sprintf('%s/KB/Carlsson_KNN_100_10_ENoise_%.3f', results_dir, sig(i)), v, 2, se_targ, 'calib_points', calib_points);        
    end
    
    %% Optimized k1 embedding  
    
    load(sprintf('%s/KB/Fit_1_10_20_10_10.mat',results_dir), 'C', 'alpha', 'theta_est', 'phi_est')
    X = arrayfun(@(t,p) DCT_proj_from_tpa(C, t, p, alpha(end-1,:)), theta_est{end-1}, phi_est{end-1}, 'un', false);
    X = vertcat(X{:});
    X = X./sqrt(sum(X.^2,2));
    [X, ~, which] = unique(X, 'rows', 'stable');
    manifold_curvature(sprintf('%s/KB/Fit_1_10_20_10_10_KNN_100_10', results_dir), X, 2, se_targ, 'calib_points', which(calib_points));
    
    %% Densely sampled image manifold
    
    spdata = load(sprintf('%s/KB/Hateren_KNN_100_10/Curvature.mat', results_dir), 'data', 'calib_points', 'ball_r');
    load(sprintf('%s/KB/HaterenAll.mat', results_dir), 'S7_unique');
    point_ids = knnsearch(S7_unique, spdata.data);
    calib_points = point_ids(calib_points);
    point_ids = unique(point_ids);
    [calib_points, where] = unique(calib_points, 'stable');
    ball_r = ball_r(where);
    manifold_curvature(sprintf('%s/KB/HaterenAll_KNN_600_10', results_dir), S7_unique, 2, se_targ, 'calib_points', calib_points, 'point_ids', point_ids);
    curvature_at_length_scale(sprintf('%s/KB/HaterenAll_KNN_600_10_DataRadius', results_dir), S7_unique, 2, ball_r, calib_points, 'point_ids', point_ids);

end