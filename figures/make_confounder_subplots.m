function make_confounder_subplots(results_dir)

    outdir = sprintf('%s/Figures/Confounder', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
    %% Non-uniform sampling
    
    nu_sparse{1} = load(sprintf('%s/Confounder/Cube_N1000_0.02/Curvature.mat', results_dir), 'S', 'dS');
    nu_sparse{2} = load(sprintf('%s/Confounder/Cube_N1000_R5/Curvature.mat', results_dir), 'S', 'dS');
    nu_sparse{3} = load(sprintf('%s/Confounder/Cube_N1000_R50/Curvature.mat', results_dir), 'S', 'dS');
    nu_sparse{4} = load(sprintf('%s/Confounder/Cube_N1000_R95/Curvature.mat', results_dir), 'S', 'dS');
    min_pct = [1 1 1 0.5];
    max_pct = [99 99 99 99];
    dS_pct  = [99.5 99 99 99];
    plot_nus_quad_error(nu_sparse, ['\color{blue}N_1=10K', '   ', '\color{green}N_2=1K'], 10000, 1000, min_pct, max_pct, dS_pct);
    paper_print(sprintf('%s/Cube_N2_1000', outdir));
        
    nu_dense{1} = load(sprintf('%s/Confounder/Cube_N10000_0.01/Curvature.mat', results_dir), 'S', 'dS');
    nu_dense{2} = load(sprintf('%s/Confounder/Cube_N10000_R5/Curvature.mat', results_dir), 'S', 'dS');
    nu_dense{3} = load(sprintf('%s/Confounder/Cube_N10000_R50/Curvature.mat', results_dir), 'S', 'dS');    
    nu_dense{4} = load(sprintf('%s/Confounder/Cube_N10000_R95/Curvature.mat', results_dir), 'S', 'dS');
    min_pct = [0.5 0.5 1 0.5];
    max_pct = [99.9 99.5 99.9 99.5];
    dS_pct  = [99 99 99 99];
    plot_nus_quad_error(nu_dense, ['\color{blue}N_1=10K', '   ', '\color{green}N_2=10K'], 10000, 10000, min_pct, max_pct, dS_pct);
    paper_print(sprintf('%s/Cube_N2_10000', outdir));

    %% Observational noise
    
    clear S;
    load(sprintf('%s/Confounder/Meta.mat', results_dir), 'S2_noise');        
    for i = 1:length(S2_noise.noise_level)
        S{i} = load(sprintf('%s/Confounder/S2_Noise_%.3f/Curvature.mat', results_dir, S2_noise.noise_level(i)), 'S', 'data', 'N_neighbors', 'dim_mfld');        
        S{i} = S{i}.S;
    end    
    plot_S2_noise_sweep(S, S2_noise.noise_level);
    paper_print(sprintf('%s/S2_Noise', outdir));    
    
    clear S2_noise
    
    %% Ambient dimension
    
    load(sprintf('%s/Confounder/Meta.mat', results_dir), 'embed');      
    dat = cell(length(embed.dims), length(embed.noise_levels));
    for i = 1:length(embed.dims)
        for j = 1:length(embed.noise_levels)        
            dat{i,j} = load(sprintf('%s/Confounder/Embed_%d_%.2f/Curvature.mat', results_dir, embed.dims(i), embed.noise_levels(j)), 'S', 'ball_r', 'data', 'N_neighbors', 'dim_mfld');                        
        end
    end    
    
    plot_Sn_noise_sweep(dat, embed.noise_levels, embed.dims);
    paper_print(sprintf('%s/Sn_Noise', outdir));    
    
    dat = cellfun(@(x) mean(x.ball_r), dat);
    plot_Sn_dim_sweep(dat, embed.noise_levels, embed.dims);
    paper_print(sprintf('%s/Sn_Dim', outdir));    
    
    clear embed

    %% Misdimension
    
    true_dim=3;    
    dat = load_misdim_data(results_dir, 'S3', true_dim);
    plot_misdim(dat, 6, true_dim, '10K Points from $\mathcal{S}^3$');
    paper_print(sprintf('%s/Misdim_S3', outdir));
        
    true_dim=4;    
    dat = load_misdim_data(results_dir, 'S2xS2', true_dim);
    plot_misdim(dat, 4, true_dim, '10K Points from $\mathcal{S}^2 \times \mathcal{S}^2$');
    paper_print(sprintf('%s/Misdim_S2xS2', outdir));

    close all;
    
end
        
function dat = load_misdim_data(results_dir, data_str, true_dim)

    dat{1} = load(sprintf('%s/Confounder/%s_Dim%d/Curvature.mat',        results_dir, data_str, true_dim-1), 'S', 'dS', 'gof', 'ball_r', 'data', 'N_neighbors', 'dim_mfld');        
    dat{2} = load(sprintf('%s/Confounder/%s_Dim%d/Curvature.mat',        results_dir, data_str, true_dim+1), 'S', 'dS', 'gof', 'ball_r', 'data', 'N_neighbors', 'dim_mfld');    
    dat{3} = load(sprintf('%s/Confounder/%s_Dim%d/Curvature.mat',        results_dir, data_str, true_dim),   'S', 'dS', 'gof', 'ball_r', 'data', 'N_neighbors', 'dim_mfld');
    
end