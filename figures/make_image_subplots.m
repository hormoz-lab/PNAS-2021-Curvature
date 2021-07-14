function make_image_subplots(results_dir)

    outdir = sprintf('%s/Figures/Image', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
    %% Schematics
    
    plot_K2_in_R3();
    paper_print(sprintf('%s/K2_in_R3', outdir));
        
    load(sprintf('%s/KB/CarlssonKB.mat', results_dir), 'Klein_R9_fcn');    
    plot_KB_patches(Klein_R9_fcn, '$\theta \in [0,\pi]$', '$\phi \in [0,2\pi]$');
    paper_print(sprintf('%s/Topology', outdir));
    
    %% Original image manifold
    
    k0 = load(sprintf('%s/KB/HaterenCarlssonMap.mat', results_dir), 'theta_map', 'phi_map');
    emp = load(sprintf('%s/KB/Hateren_KNN_100_10/Curvature.mat', results_dir), 'S', 'data', 'calib_points', 'ball_r');
    plot_S_on_2D_plane(k0.theta_map, k0.phi_map, emp.S, '$\theta$ (rad)', '$\phi$ (rad)', 'Image Patches Plotted on k^0', true);
    paper_print(sprintf('%s/Data_k0', outdir));

    %% Carlsson k0 analytical
    
    load(sprintf('%s/KB/CarlssonKB.mat', results_dir), 'Schr_fcn');
    theta_grid = linspace(0, pi, 200);
    phi_grid = linspace(0, 2*pi, 200);
    [theta_grid, phi_grid] = meshgrid(theta_grid, phi_grid);
    theta_grid = theta_grid(:);
    phi_grid = phi_grid(:);
    S = arrayfun(@(t,p) Schr_fcn(t,p), theta_grid, phi_grid);
    plot_S_on_2D_plane(theta_grid, phi_grid, S, '$\theta$ (rad)', '$\phi$ (rad)', 'Analytical k^0', false);
    paper_print(sprintf('%s/k0', outdir));    
    
    %% Carlsson k0 computational on data radius
    
    no_noise_emp = load(sprintf('%s/KB/Carlsson_KNN_100_10_Data_Radius/Curvature.mat', results_dir), 'S', 'data', 'calib_points', 'ball_r', 'dS');
    plot_S_on_2D_plane(k0.theta_map, k0.phi_map, no_noise_emp.S, '$\theta$ (rad)', '$\phi$ (rad)', {'Computed k^0 with Image', 'Patch Neighborhoods'}, true, true);
    paper_print(sprintf('%s/k0_DataRadius', outdir));    
    
    %% Carlsson k0 computational with different noise
    
    no_noise = load(sprintf('%s/KB/Carlsson_KNN_100_10_ENoise_0.000/Curvature.mat', results_dir), 'S', 'data', 'calib_points', 'ball_r', 'dS');
    plot_S_on_2D_plane(k0.theta_map, k0.phi_map, no_noise.S, '$\theta$ (rad)', '$\phi$ (rad)', 'Computed k^0', true);
    paper_print(sprintf('%s/k0_sig_0.000', outdir));    
   
    lo_noise = load(sprintf('%s/KB/Carlsson_KNN_100_10_ENoise_0.007/Curvature.mat', results_dir), 'S', 'data', 'calib_points', 'ball_r');
    plot_S_on_2D_plane(k0.theta_map, k0.phi_map, lo_noise.S, '$\theta$ (rad)', '$\phi$ (rad)', 'Computed Noisy k^0, \sigma=0.007', true);
    paper_print(sprintf('%s/k0_sig_0.007', outdir));
    
    mid_noise = load(sprintf('%s/KB/Carlsson_KNN_100_10_ENoise_0.010/Curvature.mat', results_dir), 'S', 'data', 'calib_points', 'ball_r');
    plot_S_on_2D_plane(k0.theta_map, k0.phi_map, mid_noise.S, '$\theta$ (rad)', '$\phi$ (rad)', 'Computed Noisy k^0, \sigma=0.010', true);
    paper_print(sprintf('%s/k0_sig_0.010', outdir));
    
    hi_noise = load(sprintf('%s/KB/Carlsson_KNN_100_10_ENoise_0.030/Curvature.mat', results_dir), 'S', 'data', 'calib_points', 'ball_r');
    plot_S_on_2D_plane(k0.theta_map, k0.phi_map, hi_noise.S, '$\theta$ (rad)', '$\phi$ (rad)', 'Computed Noisy k^0, \sigma=0.030', true);
    paper_print(sprintf('%s/k0_sig_0.030', outdir));

    no_noise.data = [no_noise.data zeros(size(no_noise.data,1),3)];
    del_data = sqrt(sum((no_noise.data-emp.data).^2,2));
    del_lo   = sqrt(sum((no_noise.data-lo_noise.data).^2,2));
    del_mid  = sqrt(sum((no_noise.data-mid_noise.data).^2,2));
    del_hi   = sqrt(sum((no_noise.data-hi_noise.data).^2,2));
    fprintf('Median deflection from k0 for (Data,Hi)=(%6.4f,%6.4f)\n', median(del_data), median(del_hi));
        
    %% Distribution of Euclidean distances
    
    plot_KB_distribution([del_data del_lo del_mid del_hi], {'Image Patches'; '\sigma=0.007'; '\sigma=0.010'; '\sigma=0.030'});
    paper_print(sprintf('%s/KBDistance', outdir));

    close all;
    
    %% Neighborhoods for original data and Carlsson k0 embeddings
    
    target_points = [0,        0
                     pi/2,     0
                     0,        pi/2;
                     pi/4,     pi/2;
                     pi/2,     pi/2;
                     3*pi/4,   pi/2;
                     0,        pi;
                     pi/2,     pi;
                     0,        3*pi/2;
                     pi/4,     3*pi/2;
                     pi/2,     3*pi/2;
                     3*pi/4,   3*pi/2];
                 
    plot_neighborhood_on_2D_plane(k0.theta_map, k0.phi_map, {hi_noise.data;   mid_noise.data;   lo_noise.data;   no_noise.data}, ...
                                  no_noise.calib_points,    [hi_noise.ball_r, mid_noise.ball_r, lo_noise.ball_r, no_noise.ball_r], ...
                                  target_points, '$\theta$ (rad)', '$\phi$ (rad)', 'k^0 Neighborhoods', ...
                                  {'\sigma=0.030'; '\sigma=0.010'; '\sigma=0.007'; '\sigma=0'});                                
    paper_print(sprintf('%s/KB_Neighborhoods', outdir));
                                      
    plot_neighborhood_on_2D_plane(k0.theta_map, k0.phi_map, {emp.data}, emp.calib_points, emp.ball_r, target_points, ...
                                  '$\theta$ (rad)', '$\phi$ (rad)', {'Image Patch Neighborhoods', 'Plotted on k^0'}, {'Neighborhood'});
    paper_print(sprintf('%s/Data_Neighborhoods', outdir));
                              
    %% Optimized k1 computational
    
    load(sprintf('%s/KB/Fit_1_10_10_20_10.mat', results_dir), 'theta_est', 'phi_est');
    fit = load(sprintf('%s/KB/Fit_1_10_10_20_10_KNN_100_10/Curvature.mat', results_dir));
    
    fprintf('Median deflection from k1 to Data=%6.4f\n', median(sqrt(sum((fit.data-emp.data).^2,2))));
    
    plot_S_on_2D_plane(theta_est{end-1}, phi_est{end-1}, fit.S, '$\theta$ (rad)', '$\phi$ (rad)', 'Computed k^1', true);
    paper_print(sprintf('%s/k1', outdir));
    
    %% Original data on k1 coordinates
    
    plot_S_on_2D_plane(theta_est{end-1}, phi_est{end-1}, emp.S, '$\theta$ (rad)', '$\phi$ (rad)', 'Image Patches Plotted on k^1', true);
    paper_print(sprintf('%s/Data_k1', outdir));

    %% Neighborhood for k1
    
    plot_neighborhood_on_2D_plane(k0.theta_map, k0.phi_map, {fit.data}, fit.calib_points, fit.ball_r, target_points, ...
                                  '$\theta$ (rad)', '$\phi$ (rad)', {'k^1 Neighborhoods', 'Plotted on k^0'}, {'Neighborhood'});
    paper_print(sprintf('%s/k1_Neighborhoods_on_k0', outdir));

    %% Densely sampled data on k0+k1 coordinates
    
    load(sprintf('%s/KB/SparseToFullMap.mat', results_dir));
    [idx, where] = unique(idx, 'stable');
    [idx, ind] = sort(idx);
    
    dense = load(sprintf('%s/KB/HaterenAll_KNN_600_10/Curvature.mat', results_dir), 'data', 'S', 'point_ids', 'calib_points', 'point_ids', 'ball_r');
    assert(isequal(idx, dense.point_ids));
     
    plot_S_on_2D_plane(theta_est{end-1}(where(ind)), phi_est{end-1}(where(ind)), dense.S, '$\theta$ (rad)', '$\phi$ (rad)', 'Aug. Image Patches Plotted on k^1', true);
    paper_print(sprintf('%s/Augmented_Data_k1', outdir));

    plot_neighborhood_on_2D_plane(k0.theta_map(where(ind)), k0.phi_map(where(ind)), {dense.data}, find(ismember(dense.calib_points, dense.point_ids)), dense.ball_r, target_points, ...
                                  '$\theta$ (rad)', '$\phi$ (rad)', {'Aug. Image Patch Neighborhoods', 'Plotted on k^0'}, {'Neighborhood'});
    paper_print(sprintf('%s/Augmented_Data_Neighborhoods_k0', outdir));
    
    dense_emp = load(sprintf('%s/KB/HaterenAll_KNN_600_10_DataRadius/Curvature.mat', results_dir), 'S', 'point_ids');
    assert(isequal(idx, dense_emp.point_ids));
    
    plot_S_on_2D_plane(k0.theta_map(where(ind)), k0.phi_map(where(ind)), dense_emp.S, '$\theta$ (rad)', '$\phi$ (rad)', {'Aug. Image Patches with Original', 'Neighborhoods Plotted on k^0'}, true, true);
    paper_print(sprintf('%s/Augmented_Data_Radius_k0', outdir));
        
    close all;    

end
