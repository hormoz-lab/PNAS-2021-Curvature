function make_extrinsic_subplots(results_dir)

    close all;
    
    outdir = sprintf('%s/Figures/Extrinsic', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
    %% Schematic

    plot_wavy_manifold();
    paper_print(sprintf('%s/Wavy', outdir));
    
    load(sprintf('%s/ToyData/S2/Curvature.mat', results_dir));
    
    center_id = 3;
    center_eps = 0.2;
    patch_id = rangesearch(data, data(center_id,:), center_eps);    
    local_coords = get_orthonormal_coords(data(patch_id{1},:));    
    [~, ~, ~, hessmats] = quadfit(local_coords(:,1:2), local_coords(:,3));    
    [x1, x2] = meshgrid(linspace(-center_eps, center_eps, 100), linspace(-center_eps, center_eps, 100));
    analytic = hessmats(1,1)*(x1.*x1)+hessmats(2,1)*(x1.*x2)+hessmats(3,1)*(x2.*x2);
    analytic = analytic-max(analytic(:))+max(local_coords(:,3));
    r = sqrt(x1.^2+x2.^2);
    analytic(r > center_eps) = NaN;
    
    plot_taylor_expansion(x1, x2, local_coords, analytic, center_eps);
    paper_print(sprintf('%s/Taylor', outdir));
    
    %% S2
    
    plot_manifold_with_S(data, S, [prctile(S, 1) prctile(S,99)], '10K Points from $\mathcal{S}^2$');
    paper_print(sprintf('%s/S2_Manifold', outdir));    
    plot_curvature_errorbars_2d(S, dS, 2, 1, 99, 99);
    paper_print(sprintf('%s/S2_Distribution', outdir));    
      
    %% H2
    
    load(sprintf('%s/ToyData/H2/Curvature.mat', results_dir));    
    data = data(1:length(S),:);
    plot_manifold_with_S(data, S, [prctile(S, 1) prctile(S,99)], '10K Points from $H_2^2$');
    paper_print(sprintf('%s/H2_Manifold', outdir));
    plot_curvature_with_errorbars(analytic, S, dS, 'z', data(:,3));    
    paper_print(sprintf('%s/H2_Distribution', outdir));    
       
    %% T2
    
    load(sprintf('%s/ToyData/T2/Curvature.mat', results_dir));
    plot_manifold_with_S(data, S, [prctile(S, 1) prctile(S,99)], '10K Points from $T^2$');
    paper_print(sprintf('%s/T2_Manifold', outdir));
    plot_curvature_with_errorbars(analytic, S, dS, '\theta (rad)', theta(:,1));
    paper_print(sprintf('%s/T2_Distribution', outdir));    
        
    %% S^n
    
    n = [2 3 5 7];    
    for i = 1:length(n)
        Sn{i} = load(sprintf('%s/ToyData/S%d/Curvature.mat', results_dir, n(i)), 'S');        
    end
    S = cellfun(@(x) x.S, Sn, 'UniformOutput', false);
    
    plot_Sn_distribution(n, S);
    paper_print(sprintf('%s/Sn_Distribution', outdir));
    
    close all;
    
end

