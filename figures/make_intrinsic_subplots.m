function make_intrinsic_subplots(results_dir)

    outdir = sprintf('%s/Figures/Intrinsic', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end

    load(sprintf('%s/Intrinsic.mat', results_dir), 'x', 'traces');
    
    title_str={'$z_{\infty}(x)$'; '$z_{1000}(x)$'; '$z_{49}(x)$'; '$\overline{z}_{1000}(x)$'; '$\overline{z}_{49}(x)$'; '$\tilde{z}_{49}(x;0.21)$'};    
    plot_heat_traces(x, traces, title_str);
    paper_print(sprintf('%s/HeatTrace', outdir));    
    
    load(sprintf('%s/Intrinsic.mat', results_dir), 'f', 'frac_error');
    
    estimator_names = {'Gaussian'; 'Weighted kNN (k=1000)'; 'Weighted kNN (k=100)'; 'r-Neighborhood'};
    plot_ev_fractional_error(frac_error, estimator_names, f);
    paper_print(sprintf('%s/EigenvalueError', outdir));
    
    load(sprintf('%s/Intrinsic.mat', results_dir), 'sol2', 'sol3', 'sol4');
    
    S_valid = [1.5 2.5];

    plot_poly_fit_traces(outdir, title_str, 2, x, sol2, S_valid);
    plot_poly_fit_traces(outdir, title_str, 3, x, sol3, S_valid);
    plot_poly_fit_traces(outdir, title_str, 4, x, sol4, S_valid);
    
    close all;
end

function plot_poly_fit_traces(outdir, title_str, poly_order, x, sol, S_valid)

    for i = 1:3
        sol{i}(sol{i}<S_valid(1))=NaN;
        sol{i}(sol{i}>S_valid(2))=NaN;
    end
    
    sol{4}(isnan(sol{2}))=NaN;
    sol{5}(isnan(sol{3}))=NaN;
    sol{6}(isnan(sol{3}))=NaN;
    
    fprintf('Correction fraction in D_49: %.f%%\n', sum(sol{6}>=S_valid(1) & sol{6}<=S_valid(2))/sum(~isnan(sol{3}))*100);
    
    for i = 1:length(sol)
        plot_trace_fitted_curvature(x, sol{i}, title_str{i})
        paper_print(sprintf('%s/P%d_Estimate%d', outdir, poly_order, i));
    end
    
end
