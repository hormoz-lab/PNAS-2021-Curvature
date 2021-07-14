function map_to_KB(results_dir)

    load(sprintf('%s/KB/Hateren.mat', results_dir), 'S7', 'D_knn_100');    
    S7 = S7(D_knn_100<=prctile(D_knn_100,10),:);
    [~, uniq_ind] = unique(S7, 'rows', 'stable');
    
    [C, alpha] = get_DCT_proj_matrix(false, 0, 0, 2, 1);    
    [theta_map, phi_map] = get_Klein_coords_opt(S7(uniq_ind,:), C, alpha);
    [~, uniq_angle] = unique([theta_map, phi_map], 'rows', 'stable');
    theta_map = theta_map(uniq_angle);
    phi_map = phi_map(uniq_angle);
    save(sprintf('%s/KB/HaterenCarlssonMap.mat', results_dir), 'theta_map', 'phi_map');
    
    S7_unique = S7(uniq_ind(uniq_angle),:);
    uniq_ind = uniq_ind(uniq_angle);
    save(sprintf('%s/KB/Hateren.mat', results_dir), 'S7_unique', 'uniq_ind', '-append');    
    
end