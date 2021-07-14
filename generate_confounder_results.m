function generate_confounder_results(results_dir)

    if (~exist(results_dir, 'dir'))
        mkdir(results_dir);
    end
    
    %% Non-uniform sampling
    
    rng(19193049);
    skew_sample.dim_mfld = 3;
    skew_sample.dim_amb = 8;
    skew_sample.co_dim = skew_sample.dim_amb-skew_sample.dim_mfld;
    skew_sample.se_targ = [0.02 0.01];
    skew_sample.N_dense = [1000 10000];
    skew_sample.N_uni   = 1e4;
    skew_sample.L_uni = 10;
    skew_sample.L_dense = 1;
    skew_sample.L_noise = 0.01;
    for i = 1:2        
        X = [[skew_sample.L_uni*  (rand(skew_sample.N_uni     ,skew_sample.dim_mfld)-0.5); ...
              skew_sample.L_dense*(rand(skew_sample.N_dense(i),skew_sample.dim_mfld)-0.5)], ...
              skew_sample.L_noise*randn(skew_sample.N_uni+skew_sample.N_dense(i), skew_sample.co_dim)];
        manifold_curvature(sprintf('%s/Confounder/Cube_N%d_%s', results_dir, skew_sample.N_dense(i), num2str(skew_sample.se_targ(i))), ...
                           X, skew_sample.dim_mfld, skew_sample.se_targ(i));
    end
    
    for i = 1:2
        load(sprintf('%s/Confounder/Cube_N%d_%s/Curvature.mat', results_dir, skew_sample.N_dense(i), num2str(skew_sample.se_targ(i))), ...
             'data', 'ball_r', 'which_calib');
        curvature_at_length_scale(sprintf('%s/Confounder/Cube_N%d_R5',  results_dir, skew_sample.N_dense(i)), ...
            data, skew_sample.dim_mfld, prctile(ball_r(which_calib),  5));
        curvature_at_length_scale(sprintf('%s/Confounder/Cube_N%d_R50', results_dir, skew_sample.N_dense(i)), ...
            data, skew_sample.dim_mfld, prctile(ball_r(which_calib), 50));
        curvature_at_length_scale(sprintf('%s/Confounder/Cube_N%d_R95', results_dir, skew_sample.N_dense(i)), ...
            data, skew_sample.dim_mfld, prctile(ball_r(which_calib), 95));
    end
    
    save(sprintf('%s/Confounder/Meta.mat', results_dir), 'skew_sample');
    
    %% Observational noise    
    
    N_points = 10000;    

    rng(20439);
    
    S2_noise.X = ManifoldHandler.sample_hypersphere(N_points, 2);    
    S2_noise.noise_level = exp(linspace(log(1e-3), log(3e-1), 6))';
    S2_noise.noise = randn(size(S2_noise.X));    

    for i = 1:length(S2_noise.noise_level)
        manifold_curvature(sprintf('%s/Confounder/S2_Noise_%.3f', results_dir, S2_noise.noise_level(i)), ...
                            S2_noise.X+S2_noise.noise_level(i)*S2_noise.noise, 2, 0.05);
    end
    
    save(sprintf('%s/Confounder/Meta.mat', results_dir), 'S2_noise', 'N_points', '-append');
            
    %% Ambient dimension
    
    rng(394802);
     
    embed.dims = [3,10,20,30,40,60,80,100];
    embed.noise_levels = [0.01 0.03 0.05];
        
    embed.X = ManifoldHandler.sample_hypersphere(N_points, 2);
    embed.imageplane = arrayfun(@(i) orth(rand(i,3)), embed.dims, 'un', false);
    embed.noise = arrayfun(@(i) normrnd(0,1,[N_points,i]), embed.dims, 'un', false);
    
    for i = 1:length(embed.dims)        
        for j = 1:length(embed.noise_levels)        
            fprintf('\nS2 with Dim=%d, Noise=%6.4g\n\n', embed.dims(i), embed.noise_levels(j));
            data = embed.X*embed.imageplane{i}'+embed.noise{i}*embed.noise_levels(j);
            manifold_curvature(sprintf('%s/Confounder/Embed_%d_%.2f', results_dir, embed.dims(i), embed.noise_levels(j)), data, 2, 0.05);
        end
    end
    
    save(sprintf('%s/Confounder/Meta.mat', results_dir), 'embed', '-append');
    
    %% Misdimension
    
    fprintf('Mis-dimensioning for S2xS2\n');
    
    rng(49029205);
    data = [ManifoldHandler.sample_hypersphere(N_points, 2), ManifoldHandler.sample_hypersphere(N_points, 2)];        
    misdimension(sprintf('%s/Confounder/S2xS2', results_dir), data, 4, 0.05);
                   
    fprintf('\nMis-dimensioning for S3\n\n');
                   
    rng(7777824);    
    data = ManifoldHandler.augment_dimension(ManifoldHandler.sample_hypersphere(N_points, 3), 5, 'zeropad')+0.01*randn(N_points,5);           
    misdimension(sprintf('%s/Confounder/S3', results_dir), data, 3, 0.05);
    
end

function misdimension(results_dir, data, true_dim, se_targ)

    manifold_curvature(sprintf('%s_Dim%d', results_dir, true_dim), data, true_dim, se_targ);
    manifold_curvature(sprintf('%s_Dim%d', results_dir, true_dim-1), data, true_dim-1, se_targ);
    manifold_curvature(sprintf('%s_Dim%d', results_dir, true_dim+1), data, true_dim+1, se_targ);      
    
end