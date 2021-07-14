function make_SC_subplots(data_dir, results_dir)

    outdir = sprintf('%s/Figures/SC', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    
    %% PBMC
    
    coords = csvread(sprintf('%s/PBMC/UMap.csv', data_dir), 0, 1);
    tot_counts = csvread(sprintf('%s/PBMC/TotCounts.csv', data_dir), 0, 1);
    cell_score = csvread(sprintf('%s/PBMC/CellTypeScore.csv', data_dir), 1, 1);
    [~, which] = max(cell_score,[],2);
    cell_types = {'CD4+ T cells'; 'CD8+ T cells'; 'B cells'; 'NK cells'; 'CD14+ monocytes'; 'CD16+ monocytes'; 'Dendritic cells'};
    cell_types = cell_types(which);
    target_points = [ -12   2;                    
                     -8.5  10;
                     -9.5 5.5;
                      6.5  13;
                       10  12;
                        5 -10;
                       10 -5;
                      2.5 -5;
                      4.8  2];                  
    common_subplots(results_dir, 'PBMC', coords, tot_counts, target_points, cell_types, 'PBMC (~10K cells)');         
    
    %% Gastrulation
        
    meta = readtable(sprintf('%s/Gastrulation/meta.dat', data_dir));
    PCs = readtable('Data/Gastrulation/PCA.csv');
    [is, where] = ismember(PCs.Var1, meta.cell);
    assert(all(is));
    coords = [cellfun(@str2num, meta.umapX(where)), cellfun(@str2num, meta.umapY(where))];
    tot_counts = csvread('Data/Gastrulation/TotCounts.csv', 0, 1);
    tot_counts = tot_counts(where);
    cell_types = meta.celltype(where);
    cell_types = cellfun(@(x) [upper(x(1)) lower(x(2:end))], cell_types, 'un', false);
    cell_types = strrep(cell_types, 'Erythroid', 'Erythroid ');
    cell_types = strrep(cell_types, 'Exe', 'ExE');
    cell_types = strrep(cell_types, 'Nmp', 'Neuromesodermal progenitors');
    cell_types = strrep(cell_types, 'Pgc', 'Primordial germ cells');
    cell_types = strrep(cell_types, 'Haemato', 'Haemato-');
    target_points = [9.3  7.1;
                       3  15;
                     -10 -2;
                       7 -12.5;
                       0 -15;
                     1.9  7.5;
                      12  0.4;                    
                       9 -6;
                    0.35 -4.7;
                    -5.9 -8.5;
                    -3.7  3.9;                    
                    -5.5  9];                
    common_subplots(results_dir, 'Gastrulation', coords, tot_counts, target_points, cell_types, 'Gastrulation (~120K cells)');   
    
    %% Brain
    
    coords = readtable(sprintf('%s/Brain/tsne.csv', data_dir));    
    tot_counts = csvread('Data/Brain/TotCounts.csv', 0, 1);    
    ss_cbs = readtable(sprintf('%s/Brain/cell_barcodes.csv', data_dir));
    cell_clusters = dlmread(sprintf('%s/Brain/neural_64_cluster_prediction.csv', data_dir), ',', 1, 0);
    cluster_phenotype = readtable(sprintf('%s/Brain/ClusterPhenotype.csv', data_dir));
    cell_types = cell(length(coords.Barcode),1);
    [is, where] = ismember(ss_cbs.Barcodes, coords.Barcode);
    cell_types(where(is)) = cluster_phenotype.Var2(cell_clusters(is)+1);
    coords = [coords.TSNE_1, coords.TSNE_2];
    target_points = [-6.5  33.2;
                      9.3 -30.7;
                     21.5 -20.5;
                       29 -8.5;
                      -27  5;
                    -11.5 -17.5;
                      -18  20;
                     16.5  28.5;
                      1.5  24;
                     25.5  8.5;
                     -5.5 -30;
                        9  1;
                     -4.5 -3.5];
    common_subplots(results_dir, 'Brain', coords, tot_counts, target_points, cell_types, 'Brain (~1.3M cells)');

    %% DentateGyrus
        
    coords = csvread(sprintf('%s/DentateGyrus/UMap.csv', data_dir));  
    tot_counts = csvread(sprintf('%s/DentateGyrus/TotCounts.csv', data_dir));
    cell_types = splitlines(fileread(sprintf('%s/DentateGyrus/CellType.csv', data_dir)));
    cell_types = cell_types(1:end-1);   
    cell_types = strrep(cell_types, 'Imm', 'Immature ');
    cell_types = strrep(cell_types, 'Sub', 'Subiculum');
    cell_types = strrep(cell_types, 'Radial', 'Radial ');
    cell_types = strrep(cell_types, 'Glial', 'Glial ');
    cell_types = strrep(cell_types, 'Prog', 'Progenitor');
    cell_types = strrep(cell_types, 'Nbl', 'Neuroblast');
    cell_types = strrep(cell_types, 'Astro', 'Astrocyte');
    cell_types = strrep(cell_types, '1', ' 1');
    cell_types = strrep(cell_types, '2', ' 2');
    cell_types(strcmp(cell_types, 'Radial Glia')) = {'Radial Glia 1'};
    
    target_points = [-3 0;
                     -3 6;
                      3.5 0;
                      3 6
                      9 6;
                      8 12;
                      9 18;
                      8 -4;
                      13 0;
                      15 -6];
    common_subplots(results_dir, 'DentateGyrus', coords, tot_counts, target_points, cell_types, 'Dentate Gyrus (~20K cells)');
    
end

function common_subplots(results_dir, dataset, coords, tot_counts, target_points, cell_annotations, title_str)
    
    outdir = sprintf('%s/Figures/SC/%s', results_dir, dataset);
    if (~exist(outdir, 'dir'))        
        mkdir(outdir);
    end

    plot_cell_annotation(coords, cell_annotations);
    paper_print(sprintf('%s/CellAnnotations', outdir));

    plot_value_on_UMAP(coords(:,1), coords(:,2), tot_counts, '# of Transcripts', true);
    paper_print(sprintf('%s/Transcripts', outdir));    
    
    if (isequal(dataset, 'PBMC'))
        min_pct = 0.1;
        max_pct = 99.9;
        dS_pct = 99;
    elseif (isequal(dataset, 'Gastrulation'))
        min_pct = 0.5;
        max_pct = 99.99;
        dS_pct = 99.5;
    elseif (isequal(dataset, 'Brain'))
        min_pct = 0.01;
        max_pct = 99.999;
        dS_pct = 99.99;
    elseif (isequal(dataset, 'DentateGyrus'))
        min_pct = 0.1;
        max_pct = 99.95;
        dS_pct = 99.5;
    else
       min_pct = 0;
       max_pct = 100;
       dS_pct = 100;
    end
    
    load(sprintf('%s/SC/%s_Meta.mat', results_dir, dataset), 'knn', 'ss_id', 'ss_transcript_id');
    
    R25 = load(sprintf('%s/SC/%s_R25/Curvature.mat', results_dir, dataset), 'S', 'ball_r');
    R50 = load(sprintf('%s/SC/%s_R50/Curvature.mat', results_dir, dataset), 'S', 'ball_r');
    R75 = load(sprintf('%s/SC/%s_R75/Curvature.mat', results_dir, dataset), 'S', 'ball_r');
       
    R25.S = nanmean(R25.S(knn),2);
    R50.S = nanmean(R50.S(knn),2);
    R75.S = nanmean(R75.S(knn),2);
    
    S_min = prctile([R25.S; R50.S; R75.S],  1);
    S_max = prctile([R25.S; R50.S; R75.S], 99);
    
    plot_value_on_UMAP(coords(:,1), coords(:,2), R25.S, 'r_{25}', [S_min, S_max]);
    paper_print(sprintf('%s/S_R25_UMap', outdir));    
    plot_value_on_UMAP(coords(:,1), coords(:,2), R50.S, 'r_{50}', [S_min, S_max]);
    paper_print(sprintf('%s/S_R50_UMap', outdir));    
    plot_value_on_UMAP(coords(:,1), coords(:,2), R75.S, 'r_{75}', [S_min, S_max]);
    paper_print(sprintf('%s/S_R75_UMap', outdir));
    
    close all;
    
    pos = load(sprintf('%s/SC/%s_DimPos/Curvature.mat', results_dir, dataset), 'S', 'dim_mfld');
    plot_value_on_UMAP(coords(:,1), coords(:,2), nanmean(pos.S(knn),2), sprintf('Scalar Curvature, d=%d', pos.dim_mfld), true);
    paper_print(sprintf('%s/SPos_UMap', outdir));  
    
    neg = load(sprintf('%s/SC/%s_DimNeg/Curvature.mat', results_dir, dataset), 'S', 'dim_mfld');
    plot_value_on_UMAP(coords(:,1), coords(:,2), nanmean(neg.S(knn),2), sprintf('Scalar Curvature, d=%d', neg.dim_mfld), true);
    paper_print(sprintf('%s/SNeg_UMap', outdir));    
    
    clear pos neg;
    close all;
    
    dat = load(sprintf('%s/SC/%s/Curvature.mat', results_dir, dataset), ...
               'S', 'dS', 'gof', 'N_neighbors', 'data', 'which_calib', 'calib_points', 'ball_r');
    assert(size(dat.data,1) == length(coords));

    plot_cell_type_curvatures(cell_annotations, dat.S, dat.N_neighbors, 'Global');
    paper_print(sprintf('%s/Ttest_Cell_PW', outdir));    
    
    plot_cell_type_boxplots(cell_annotations, dat.S, dataset);
    paper_print(sprintf('%s/Boxplot', outdir));    
    
    close all;
    
    plot_value_on_UMAP(coords(:,1), coords(:,2), nanmean(dat.S(knn),2), title_str, true); 
    paper_print(sprintf('%s/S_UMap', outdir));    
    plot_curvature_errorbars_2d(dat.S, dat.dS, 0, min_pct, max_pct, dS_pct);
    paper_print(sprintf('%s/S_Err', outdir));
    
    plot_value_on_UMAP(coords(:,1), coords(:,2), dat.gof, 'GOF p-value', false);
    paper_print(sprintf('%s/GOF_UMap', outdir));    
    plot_stat_histogram(dat.gof, 'GOF p-value', false);
    paper_print(sprintf('%s/GOF_Histogram', outdir));
    
    plot_value_on_UMAP(coords(:,1), coords(:,2), dat.N_neighbors, '# of Points in Neighborhood', true); 
    paper_print(sprintf('%s/Neighbors_UMap', outdir));    
    plot_stat_histogram(dat.N_neighbors, '# of Points in Neighborhood', true, nanmedian(dat.N_neighbors));
    paper_print(sprintf('%s/Neighbors_Histogram', outdir));    
    plot_neighborhood_on_UMAP(coords(:,1), coords(:,2), dat.data, dat.calib_points, dat.ball_r, target_points);
    paper_print(sprintf('%s/Neighborhood', outdir));
                
    plot_value_on_UMAP(coords(:,1), coords(:,2), dat.ball_r(dat.which_calib), 'Neighborhood Size, r', true); 
    paper_print(sprintf('%s/Radius_UMap', outdir));    
    plot_stat_histogram(dat.ball_r(dat.which_calib), 'Neighborhood Size, r', true, [R25.ball_r, R50.ball_r, R75.ball_r]);
    paper_print(sprintf('%s/Radius_Histogram', outdir));
    
    clear R25 R50 R75;
    close all;
        
    SS2 = load(sprintf('%s/SC/%s_SS2/Curvature.mat', results_dir, dataset), 'S');
    SS4 = load(sprintf('%s/SC/%s_SS4/Curvature.mat', results_dir, dataset), 'S');
    
    plot_downsample_correlation(dat.S(ss_id{1}), SS2.S, 'Original', 'Downsampled', {'Cells Downsampled'; '2x (f=50%)'}, 'Pearson', false, [1, 99]);    
    paper_print(sprintf('%s/CellDownsampleCorrelation_SSFact2', outdir));

    plot_downsample_correlation(dat.S(ss_id{2}), SS4.S, 'Original', 'Downsampled', {'Cells Downsampled'; '4x (f=75%)'}, 'Pearson', false, [1, 99]);    
    paper_print(sprintf('%s/CellDownsampleCorrelation_SSFact4', outdir));
        
    [~, mask] = ismember(knn(ss_id{1},:), ss_id{1});
    SS2 = arrayfun(@(i) nanmean(SS2.S(nonzeros(mask(i,:)))), [1:length(SS2.S)]');
    
    [~, mask] = ismember(knn(ss_id{2},:), ss_id{2});
    SS4 = arrayfun(@(i) nanmean(SS4.S(nonzeros(mask(i,:)))), [1:length(SS4.S)]');
    
    S_min = prctile([SS2; SS4],  1);
    S_max = prctile([SS2; SS4], 99);
    
    plot_value_on_UMAP(coords(ss_id{1},1), coords(ss_id{1},2), SS2, 'Cells Downsampled 2x (f=50%)', [S_min, S_max]);
    paper_print(sprintf('%s/S_SSFact2_UMap', outdir));    
    
    plot_value_on_UMAP(coords(ss_id{2},1), coords(ss_id{2},2), SS4, 'Cells Downsampled 4x (f=75%)', [S_min, S_max]);
    paper_print(sprintf('%s/S_SSFact4_UMap', outdir));        
    
    if (isequal(dataset, 'PBMC'))
        
        SS2 = load(sprintf('%s/SC/%s_SSTranscript2/Curvature.mat', results_dir, dataset), 'S');
        SS4 = load(sprintf('%s/SC/%s_SSTranscript4/Curvature.mat', results_dir, dataset), 'S');
        
        plot_downsample_correlation(dat.S(ss_transcript_id{1}), SS2.S, 'Original', 'Downsampled', {'Transcripts Downsampled'; '2x (f=50%)'}, 'Pearson', false, [1, 99]);        
        paper_print(sprintf('%s/TranscriptDownsampleCorrelation_SSFact2', outdir));

        plot_downsample_correlation(dat.S(ss_transcript_id{2}), SS4.S, 'Original', 'Downsampled', {'Transcripts Downsampled'; '4x (f=75%)'}, 'Pearson', false, [1, 99]);        
        paper_print(sprintf('%s/TranscriptDownsampleCorrelation_SSFact4', outdir));
        
        [~, mask] = ismember(knn(ss_transcript_id{1},:), ss_transcript_id{1});
        SS2 = arrayfun(@(i) nanmean(SS2.S(nonzeros(mask(i,:)))), [1:length(SS2.S)]');
    
        [~, mask] = ismember(knn(ss_transcript_id{2},:), ss_transcript_id{2});
        SS4 = arrayfun(@(i) nanmean(SS4.S(nonzeros(mask(i,:)))), [1:length(SS4.S)]');
        
        S_min = prctile([SS2; SS4],  1);
        S_max = prctile([SS2; SS4], 99);
    
        plot_value_on_UMAP(coords(ss_transcript_id{1},1), coords(ss_transcript_id{1},2), SS2, {'Transcripts Downsampled 2x', '(f=50%)'}, [S_min, S_max]);
        paper_print(sprintf('%s/S_SSTranscriptFact2_UMap', outdir));    
    
        plot_value_on_UMAP(coords(ss_transcript_id{2},1), coords(ss_transcript_id{2},2), SS4, {'Transcripts Downsampled 4x', '(f=75%)'}, [S_min, S_max]);
        paper_print(sprintf('%s/S_SSTranscriptFact4_UMap', outdir));        
       
    end
        
    clear SS2 SS4;
    close all;
    
    load(sprintf('%s/SC/%s_ErrorAnalysis.mat', results_dir, dataset));
    
    plot_kNN_stat(rho_pos, rho_neg, which_K, 'k', 'Correlation with kNN', 'Spatial Correlation in Scalar Curvature');
    paper_print(sprintf('%s/SpatialCorrelation', outdir));   
    
    plot_kNN_stat(frac_in_eb_pos*100, frac_in_eb_neg*100, which_K, 'k', '% Points w/ kNN in 95% CI', 'Spatial Precision of Errorbars');
    paper_print(sprintf('%s/SpatialErrorbars', outdir));   
    
    plot_SS_stat(frac_in_SS_pos*100, frac_in_SS_neg*100, [50 75], 'Downsample Factor, f (%)', '% Points in 95% CI', {'Sensitivity to'; 'Cell Downsampling'});
    paper_print(sprintf('%s/DownsamplingCell', outdir));       
    
    if (isequal(dataset, 'PBMC'))
        
        plot_SS_stat(frac_in_SS_transcript_pos*100, frac_in_SS_transcript_neg*100, [50 75], 'Downsample Factor, f (%)', '% Points in 95% CI', {'Sensitivity to'; 'Transcript Downsampling'});
        paper_print(sprintf('%s/DownsamplingTranscript', outdir));   
        
    end

    close all;

end