function plot_Sn_dim_sweep(dat, noise_levels, dims)
    
    assert(size(dat,1) == length(dims));
    assert(size(dat,2) == length(noise_levels));
    assert(length(noise_levels) == 3);
    
    fig_width = 3.5;
    fig_height = 4;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
    
    sp1 = subplot(1,1,1);
    hold on;
    for i = 1:3        
        plot(dims, dat(:,i), 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 15);
    end
    hold off;
    xticks(dims);
    xticklabels(cellstr(num2str(dims')));

    box off;            
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');        
    xlabel('Ambient Dimension, n', 'FontSize', 6);        
    ylabel('Average Ball Radius, r', 'FontSize', 6);
    
    lgd = legend(cellstr(num2str(noise_levels', '%#-.2f')), 'FontSize', 5, 'Location', 'SouthOutside', 'NumColumns', 2);
    title(lgd, '\sigma', 'FontSize', 6, 'FontWeight', 'normal');
    legend('boxoff');
    
    left_margin = 0.7;
    right_margin = 0.2;    
    top_margin = 0.1;
    bottom_margin = 1.6;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
    
    lgd.Units = 'centimeters';
    lgd_height = 1.0;
    lgd.Position = [0 0 fig_width lgd_height];
        
end
