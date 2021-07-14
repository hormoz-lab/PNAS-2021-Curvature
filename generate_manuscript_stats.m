function generate_manuscript_stats(data_dir, results_dir)
    
    fid = fopen(sprintf('%s/ManuscriptStats.txt', results_dir), 'wt');
    sig_level = 0.05;
    
    %% Extrinsic
    
    fprintf(fid, 'EXTRINSIC\n\n');
    
    load(sprintf('%s/ToyData/S%d/Curvature.mat', results_dir, 2));
    analytic = ManifoldHandler.hypersphere_curvature(2);    
    fprintf(fid, 'S2: Neighborhood Population=%6.2f+/-%6.2f, Mean Error=%6.2f\n', nanmean(N_neighbors), nanstd(N_neighbors), nanmean(abs(S-analytic)));    
    fprintf(fid, '    Accuracy=%6.2f, GOF(<=%.2f)=%6.2f, sigma_h=%6.4f\n\n', pct_accurate(S, dS, analytic), sig_level, pct_fit(gof, sig_level), se_targ);
    
    load(sprintf('%s/ToyData/H2/Curvature.mat', results_dir));    
    fprintf(fid, 'H2: Accuracy=%6.2f, GOF(<=%.2f)=%6.2f, sigma_h=%6.4f\n\n', pct_accurate(S, dS, analytic), sig_level, pct_fit(gof, sig_level), se_targ);
    
    load(sprintf('%s/ToyData/T2/Curvature.mat', results_dir));    
    fprintf(fid, 'T2: Accuracy=%6.2f, GOF(<=%.2f)=%6.2f, sigma_h=%6.4f\n\n', pct_accurate(S, dS, analytic), sig_level, pct_fit(gof, sig_level), se_targ);
    
    n = [3 5 7];
    
    for i = 1:length(n)
        load(sprintf('%s/ToyData/S%d/Curvature.mat', results_dir, n(i)));
        analytic = ManifoldHandler.hypersphere_curvature(n(i));
        fprintf(fid, 'S%d: Accuracy=%6.2f, GOF(<=%.2f)=%6.2f, sigma_h=%6.4f\n\n', n(i), pct_accurate(S, dS, analytic), sig_level, pct_fit(gof, sig_level), se_targ);
    end
    
    %% Confounder
    
    fprintf(fid, 'CONFOUNDER\n\n');
    
    load(sprintf('%s/Confounder/Meta.mat', results_dir), 'skew_sample');
    
    for i = 1:length(skew_sample.se_targ)
        load(sprintf('%s/Confounder/Cube_N%d_%s/Curvature.mat', results_dir, skew_sample.N_dense(i), num2str(skew_sample.se_targ(i))), 'gof');
        fprintf(fid, 'Non-Uniform Sampling: N_dense=%6d, GOF(<=%.2f)=%6.2f, sigma_h=%6.4f\n\n', ...
            skew_sample.N_dense(i), sig_level, pct_fit(gof, sig_level), skew_sample.se_targ(i));
    end
    
    %% Image Patch and KB
    
    fprintf(fid, 'IMAGES AND KB\n\n');
    
    k0 = load(sprintf('%s/KB/HaterenCarlssonMap.mat', results_dir), 'theta_map', 'phi_map');
    load(sprintf('%s/KB/CarlssonKB.mat', results_dir), 'Schr_fcn');
    no_noise = load(sprintf('%s/KB/Carlsson_KNN_100_10_ENoise_0.000/Curvature.mat', results_dir), 'S', 'dS', 'gof', 'se_targ', 'data');
    analytic = arrayfun(@(t,p) Schr_fcn(t,p), k0.theta_map, k0.phi_map);
    fprintf(fid, 'k0: Accuracy=%6.2f, GOF(<=%.2f)=%6.2f, sigma_h=%6.4f\n\n', ...
        pct_accurate(no_noise.S, no_noise.dS, analytic), sig_level, pct_fit(gof, sig_level), no_noise.se_targ);
    
    no_noise.data = [no_noise.data zeros(size(no_noise.data,1),3)];
    hi_noise = load(sprintf('%s/KB/Carlsson_KNN_100_10_ENoise_0.030/Curvature.mat', results_dir), 'data');
    emp = load(sprintf('%s/KB/Hateren_KNN_100_10/Curvature.mat', results_dir), 'data');
    fit = load(sprintf('%s/KB/Fit_1_10_10_20_10_KNN_100_10/Curvature.mat', results_dir), 'data');
    
    fprintf(fid, 'Median deflection (k0, Data      )=%6.4f\n', med_dist(no_noise.data, emp.data));
    fprintf(fid, 'Median deflection (k0, High Noise)=%6.4f\n', med_dist(no_noise.data, hi_noise.data));
    fprintf(fid, 'Median deflection (k1, Data      )=%6.4f\n\n', med_dist(fit.data, emp.data));
    
    %% Single-Cell
    
    fprintf(fid, 'SC-DATA\n\n');
    
    summarize_sc(results_dir', 'PBMC', fid, sig_level, 250, 0.75);
    summarize_sc(results_dir', 'Gastrulation', fid, sig_level, 250, 0.75);
    summarize_sc(results_dir', 'Brain', fid, sig_level, 250, 0.75);
    
    load(sprintf('%s/SC/PBMC/LinearModel.mat', results_dir), 'sig_genes');
    fprintf(fid, 'PBMC: Significant Genes\n');
    fprintf(fid, '%s\n', sig_genes{:});
    fprintf(fid, '\n');
    
    load(sprintf('%s/SC/DentateGyrus/LinearModel.mat', results_dir), 'sig_genes');
    fprintf(fid, 'Dentate Gyrus: Significant Genes\n');
    fprintf(fid, '%s\n', sig_genes{:});
    fprintf(fid, '\n');
    
    load(sprintf('%s/SC/DentateGyrus/Curvature.mat', results_dir), 'S');
    
    pca_speed = csvread(sprintf('%s/DentateGyrus/PCA_Speed.csv', data_dir));
    [R, p] = corr(S, pca_speed, 'Rows', 'complete');
    fprintf(fid, 'Dentate Gyrus: Speed Correlation: rho=%6.4f, p=%8.6f\n', R, p);
    
    pca_diver = csvread(sprintf('%s/DentateGyrus/PCA_Divergence.csv', data_dir));
    [R, p] = corr(S, pca_diver, 'Rows', 'complete');
    fprintf(fid, 'Dentate Gyrus: Divergence Correlation: rho=%6.4f, p=%8.6f\n', R, p);   
    
    
end

function summarize_sc(results_dir, which, fid, sig_level, k, ss_thresh)

    dim = load(sprintf('%s/SC/%s/Curvature.mat', results_dir, which), 'data', 'dim_mfld', 'se_targ', 'gof');    
    neg = load(sprintf('%s/SC/%s_DimNeg/Curvature.mat', results_dir, which), 'se_targ', 'gof', 'dim_mfld');
    pos = load(sprintf('%s/SC/%s_DimPos/Curvature.mat', results_dir, which), 'se_targ', 'gof', 'dim_mfld');    
    
    err = load(sprintf('%s/SC/%s_ErrorAnalysis.mat', results_dir, which));
    
    fprintf(fid, '%s\n', which);
    fprintf(fid, 'N=%d, n=%d, L=%6.4f\n', size(dim.data, 1), size(dim.data,2), get_length_scale(dim.data, dim.dim_mfld));
    fprintf(fid, 'd=%d, sigma_h=%6.4f, GOF(<=%.2f)=%6.2f\n', neg.dim_mfld, neg.se_targ, sig_level, pct_fit(neg.gof, sig_level));
    fprintf(fid, 'd=%d, sigma_h=%6.4f, GOF(<=%.2f)=%6.2f\n', dim.dim_mfld, dim.se_targ, sig_level, pct_fit(dim.gof, sig_level));
    fprintf(fid, 'd=%d, sigma_h=%6.4f, GOF(<=%.2f)=%6.2f\n', pos.dim_mfld, pos.se_targ, sig_level, pct_fit(pos.gof, sig_level));
        
    err.rho_pos        = err.rho_pos       (err.which_K==k);    
    err.p_pos          = err.p_pos         (err.which_K==k);
    err.frac_in_eb_pos = err.frac_in_eb_pos(err.which_K==k);
    fprintf(fid, 'Spatial Correlation with k=%2d: rho=%6.4f, p=%8.6f\n', k, err.rho_pos, err.p_pos);
    fprintf(fid, 'Spatial Errorbar with k=%2d: Accuracy=%6.4f%%, p=%8.6f\n', k, err.frac_in_eb_pos*100, perm_test_p_val(err.frac_in_eb_pos, err.frac_in_eb_neg));
    
    ind = find(err.frac_in_SS_pos>=ss_thresh, 1, 'last');
    err.frac_in_SS_pos = err.frac_in_SS_pos(ind);
    err.frac_in_SS_neg = err.frac_in_SS_neg(ind,:);
    
    fprintf(fid, 'Cell Subsample: f=%3.2f, Accuracy=%6.4f%%, p=%8.6f\n', 1-2^-ind, err.frac_in_SS_pos*100, perm_test_p_val(err.frac_in_SS_pos, err.frac_in_SS_neg));
    
    if (strcmp(which, 'PBMC'))
        fprintf(fid, 'Transcript Subsample: f=0.50, Accuracy=%6.4f%%, p=%8.6f\n', err.frac_in_SS_transcript_pos(1)*100, perm_test_p_val(err.frac_in_SS_transcript_pos(1), err.frac_in_SS_transcript_neg(1,:)));
        fprintf(fid, '                      f=0.75, Accuracy=%6.4f%%, p=%8.6f\n', err.frac_in_SS_transcript_pos(2)*100, perm_test_p_val(err.frac_in_SS_transcript_pos(2), err.frac_in_SS_transcript_neg(2,:)));
    end
    fprintf(fid, '\n');
    
end

function pct = pct_accurate(S, dS, analytic)
    pct = sum((S+2*dS)>analytic & (S-2*dS)<analytic)/sum(~isnan(S))*100;
end

function pct = pct_fit(gof, sig_level)
    pct = sum(gof<=sig_level)/sum(~isnan(gof))*100;
end

function d = med_dist(X1, X2)
    d = median(sqrt(sum((X1-X2).^2,2)));
end

function L = get_length_scale(X, dim_mfld)
    L = 3*std(X(:,dim_mfld));
end

function p = perm_test_p_val(test, ref)
    p = sum(ref>test)/length(ref);
end