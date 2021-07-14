function generate_SC_results(data_path, results_dir)

    %% PBMC

    data = csvread(sprintf('%s/PBMC/PCA.csv', data_path), 0, 1);
    [dim_amb, dim_mfld, L] = set_params(data);
    fprintf('PBMC: %d %d %.3f\n', dim_amb, dim_mfld, L);
    data = data(:,1:dim_amb);
    
    manifold_curvature(sprintf('%s/SC/PBMC',        results_dir), data, dim_mfld  , 0.041);
    manifold_curvature(sprintf('%s/SC/PBMC_DimPos', results_dir), data, dim_mfld+1, 0.045);
    manifold_curvature(sprintf('%s/SC/PBMC_DimNeg', results_dir), data, dim_mfld-1, 0.031);
    
    load(sprintf('%s/SC/PBMC/Curvature.mat', results_dir), 'data', 'dim_mfld', 'ball_r', 'which_calib');
    
    curvature_at_length_scale(sprintf('%s/SC/PBMC_R25', results_dir), data, dim_mfld, prctile(ball_r(which_calib), 25));
    curvature_at_length_scale(sprintf('%s/SC/PBMC_R50', results_dir), data, dim_mfld, prctile(ball_r(which_calib), 50));
    curvature_at_length_scale(sprintf('%s/SC/PBMC_R75', results_dir), data, dim_mfld, prctile(ball_r(which_calib), 75));

    knn = knnsearch_block(data, 250);
    save(sprintf('%s/SC/PBMC_Meta.mat', results_dir), 'knn');

    ss_id = get_subsample_index(size(data,1));    
    save(sprintf('%s/SC/PBMC_Meta.mat', results_dir), 'ss_id', '-append');
    curvature_at_length_scale(sprintf('%s/SC/PBMC_SS2', results_dir), data(ss_id{1},:), dim_mfld, ball_r(which_calib(ss_id{1})), [1:length(ss_id{1})]');
    curvature_at_length_scale(sprintf('%s/SC/PBMC_SS4', results_dir), data(ss_id{2},:), dim_mfld, ball_r(which_calib(ss_id{2})), [1:length(ss_id{2})]');    

    PCs_orig = readtable(sprintf('%s/PBMC/PCA.csv', data_path));
    PCs_SS2 = readtable(sprintf('%s/PBMC/SSTranscripts_Fact1/PCA.csv', data_path));
    PCs_SS4 = readtable(sprintf('%s/PBMC/SSTranscripts_Fact2/PCA.csv', data_path));
    [is1, where1] = ismember(PCs_SS2.Var1, PCs_orig.Var1);
    [is2, where2] = ismember(PCs_SS4.Var1, PCs_orig.Var1);
    ss_transcript_id{1} = where1(is1);
    ss_transcript_id{2} = where2(is2);
    save(sprintf('%s/SC/PBMC_Meta.mat', results_dir), 'ss_transcript_id', '-append');

    curvature_at_length_scale(sprintf('%s/SC/PBMC_SSTranscript2', results_dir), table2array(PCs_SS2(is1,2:size(data,2)+1)), ...
        dim_mfld, ball_r(which_calib(ss_transcript_id{1})), [1:length(ss_transcript_id{1})]');
    curvature_at_length_scale(sprintf('%s/SC/PBMC_SSTranscript4', results_dir), table2array(PCs_SS4(is2,2:size(data,2)+1)), ...
        dim_mfld, ball_r(which_calib(ss_transcript_id{2})), [1:length(ss_transcript_id{2})]');
    
    %% Gastrulation

    data = csvread(sprintf('%s/Gastrulation/PCA.csv', data_path), 1, 1);    
    [dim_amb, dim_mfld, L] = set_params(data);
    fprintf('Gastrulation: %d %d %.3f\n', dim_amb, dim_mfld, L);
    data = data(:,1:dim_amb);
    
    manifold_curvature(sprintf('%s/SC/Gastrulation',        results_dir), data, dim_mfld  , 0.044);
    manifold_curvature(sprintf('%s/SC/Gastrulation_DimPos', results_dir), data, dim_mfld+1, 0.053);
    manifold_curvature(sprintf('%s/SC/Gastrulation_DimNeg', results_dir), data, dim_mfld-1, 0.036);
    
    load(sprintf('%s/SC/Gastrulation/Curvature.mat', results_dir), 'data', 'dim_mfld', 'ball_r', 'which_calib');
    
    curvature_at_length_scale(sprintf('%s/SC/Gastrulation_R25', results_dir), data, dim_mfld, prctile(ball_r(which_calib), 25));
    curvature_at_length_scale(sprintf('%s/SC/Gastrulation_R50', results_dir), data, dim_mfld, prctile(ball_r(which_calib), 50));
    curvature_at_length_scale(sprintf('%s/SC/Gastrulation_R75', results_dir), data, dim_mfld, prctile(ball_r(which_calib), 75));

    knn = knnsearch_block(data, 250);
    save(sprintf('%s/SC/Gastrulation_Meta.mat', results_dir), 'knn');

    ss_id = get_subsample_index(size(data,1));    
    save(sprintf('%s/SC/Gastrulation_Meta.mat', results_dir), 'ss_id', '-append');
    curvature_at_length_scale(sprintf('%s/SC/Gastrulation_SS2', results_dir), data(ss_id{1},:), dim_mfld, ball_r(which_calib(ss_id{1})), [1:length(ss_id{1})]');
    curvature_at_length_scale(sprintf('%s/SC/Gastrulation_SS4', results_dir), data(ss_id{2},:), dim_mfld, ball_r(which_calib(ss_id{2})), [1:length(ss_id{2})]');

    %% Brain

    data = csvread(sprintf('%s/Brain/PCA.csv', data_path), 1, 1);    
    [dim_amb, dim_mfld, L] = set_params(data);    
    fprintf('Brain: %d %d %.3f\n', dim_amb, dim_mfld, L);
    data = data(:,1:dim_amb);
    
    manifold_curvature(sprintf('%s/SC/Brain',        results_dir), data, dim_mfld  , 0.050);
    manifold_curvature(sprintf('%s/SC/Brain_DimPos', results_dir), data, dim_mfld+1, 0.055);
    manifold_curvature(sprintf('%s/SC/Brain_DimNeg', results_dir), data, dim_mfld-1, 0.034);
    
    load(sprintf('%s/SC/Brain/Curvature.mat', results_dir), 'data', 'dim_mfld', 'ball_r', 'which_calib');
    
    curvature_at_length_scale(sprintf('%s/SC/Brain_R25', results_dir), data, dim_mfld, prctile(ball_r(which_calib), 25));
    curvature_at_length_scale(sprintf('%s/SC/Brain_R50', results_dir), data, dim_mfld, prctile(ball_r(which_calib), 50));
    curvature_at_length_scale(sprintf('%s/SC/Brain_R75', results_dir), data, dim_mfld, prctile(ball_r(which_calib), 75));

    knn = knnsearch_block(data, 250);
    save(sprintf('%s/SC/Brain_Meta.mat', results_dir), 'knn', '-v7.3', '-nocompression');
    
    ss_id = get_subsample_index(size(data,1));    
    save(sprintf('%s/SC/Brain_Meta.mat', results_dir), 'ss_id', '-append');
    curvature_at_length_scale(sprintf('%s/SC/Brain_SS2', results_dir), data(ss_id{1},:), dim_mfld, ball_r(which_calib(ss_id{1})), [1:length(ss_id{1})]');
    curvature_at_length_scale(sprintf('%s/SC/Brain_SS4', results_dir), data(ss_id{2},:), dim_mfld, ball_r(which_calib(ss_id{2})), [1:length(ss_id{2})]');    

    %% DentateGyrus
    
    data = csvread(sprintf('%s/DentateGyrus/PCA.csv', data_path), 0, 0);    
    [dim_amb, dim_mfld, L] = set_params(data);
    fprintf('DentateGyrus: %d %d %.3f\n', dim_amb, dim_mfld, L);
    data = data(:,1:dim_amb);
    
    manifold_curvature(sprintf('%s/SC/DentateGyrus',        results_dir), data, dim_mfld  , 0.053);
    manifold_curvature(sprintf('%s/SC/DentateGyrus_DimNeg', results_dir), data, dim_mfld+1, 0.049);
    manifold_curvature(sprintf('%s/SC/DentateGyrus_DimPos', results_dir), data, dim_mfld+2, 0.061);
    
    load(sprintf('%s/SC/DentateGyrus/Curvature.mat', results_dir), 'data', 'dim_mfld', 'ball_r', 'which_calib');
    
    curvature_at_length_scale(sprintf('%s/SC/DentateGyrus_R25', results_dir), data, dim_mfld, prctile(ball_r(which_calib), 25));
    curvature_at_length_scale(sprintf('%s/SC/DentateGyrus_R50', results_dir), data, dim_mfld, prctile(ball_r(which_calib), 50));
    curvature_at_length_scale(sprintf('%s/SC/DentateGyrus_R75', results_dir), data, dim_mfld, prctile(ball_r(which_calib), 75));

    knn = knnsearch_block(data, 250);
    save(sprintf('%s/SC/DentateGyrus_Meta.mat', results_dir), 'knn');

    ss_id = get_subsample_index(size(data,1));    
    save(sprintf('%s/SC/DentateGyrus_Meta.mat', results_dir), 'ss_id', '-append');
    curvature_at_length_scale(sprintf('%s/SC/DentateGyrus_SS2', results_dir), data(ss_id{1},:), dim_mfld, ball_r(which_calib(ss_id{1})), [1:length(ss_id{1})]');
    curvature_at_length_scale(sprintf('%s/SC/DentateGyrus_SS4', results_dir), data(ss_id{2},:), dim_mfld, ball_r(which_calib(ss_id{2})), [1:length(ss_id{2})]');
    
end

function [dim_amb, dim_mfld, L] = set_params(data)

    vars = var(data, [], 1);
    dim_amb  = find(cumsum(vars)/sum(vars) <= 0.8, 1, 'last');
    dim_mfld = find(cumsum(vars)/sum(vars) <= 0.64, 1, 'last');
    L = 3*sqrt(vars(dim_mfld));    
end

function ss_id = get_subsample_index(N)

    rng(349343);
    ss_id{1} = datasample([1:N]'  , round(N/2)               , 'Replace', false);
    ss_id{2} = datasample(ss_id{1}, round(length(ss_id{1})/2), 'Replace', false);
    ss_id{3} = datasample(ss_id{2}, round(length(ss_id{2})/2), 'Replace', false);    
    
end