function plot_cell_type_boxplots(cell_type, S, title_str)

    unannotated = cellfun(@isempty, cell_type);
    
    cell_type(unannotated,:) = [];
    S(unannotated,:) = [];    

    [cell_type, ~, which_type] = unique(cell_type, 'stable');    
    [~, cell_order] = sort(accumarray(which_type, 1), 'descend');
    cell_type = cell_type(cell_order);
    [~, which_type] = ismember(which_type, cell_order);
    N_type = length(cell_type);
    
    q1 = arrayfun(@(i) quantile(S(which_type==i), 0.25), [1:N_type]');
    q3 = arrayfun(@(i) quantile(S(which_type==i), 0.75), [1:N_type]');
    
    pad = 0.0;
    y_min = min(q1-1.5*(q3-q1))-pad;
    y_max = max(q3+1.5*(q3-q1))+pad;
    
    left_margin = 0.7;
    right_margin = 0.1;
    top_margin = 0.3;
    bottom_margin = 2.75;
    
    rel_width = 0.25*N_type;
    rel_height = 3;
    
    fig_width = left_margin+right_margin+rel_width;
    fig_height = bottom_margin+top_margin+rel_height;
        
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
        
    sp1 = subplot(1,1,1);
    
    boxplot(S, which_type,'PlotStyle','compact', 'OutlierSize', 0.1, 'Labels', cell_type, 'Symbol', '');
    ylim([y_min y_max]);
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca,'Xtick', [1:N_type], 'XTickLabel',cell_type,'FontSize',5)
    xtickangle(90);
    ylabel('Scalar Curvature, S', 'FontSize', 6);
    title(title_str, 'FontSize', 6, 'FontWeight', 'Normal');
        
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin bottom_margin rel_width rel_height]);

end