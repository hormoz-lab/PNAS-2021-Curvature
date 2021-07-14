function fit_curvature_to_gene_exp(data_dir, results_dir, dataset)

    genes = splitlines(fileread(sprintf('%s/%s/PCA_Genes.csv', data_dir, dataset)));
    genes = genes(1:end-1);

    gene_exp = csvread(sprintf('%s/%s/NormalizedCounts.csv', data_dir, dataset));    
    
    load(sprintf('%s/SC/%s/Curvature.mat', results_dir, dataset), 'S');
    
    mdl = fitlm(gene_exp, S, 'VarNames', [genes; 'S']);    
    fdr_rate = 0.05;    
    p_crit = benjamini_hochberg(mdl.Coefficients.pValue(2:end), fdr_rate);
    
    sig_genes = mdl.CoefficientNames(1+find(mdl.Coefficients.pValue(2:end)<=p_crit));
    
    save(sprintf('%s/SC/%s/LinearModel.mat', results_dir, dataset), 'mdl', 'fdr_rate', 'p_crit', 'sig_genes');
    
end