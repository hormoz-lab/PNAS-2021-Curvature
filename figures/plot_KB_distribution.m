function plot_KB_distribution(del, series_names)

    assert(size(del,2) == length(series_names));

    fig_width = 3.5;
    fig_height = 3.4;
        
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
       
    sp1 = subplot(1,1,1);
    hold on;
    
    for i = 1:size(del,2)
        histogram(del(:,i),[0:0.01:1.0],'Normalization','pdf', 'DisplayName', series_names{i});
    end    
    
    hold off;
    axis tight;
    xlim([0 0.5]);
    
    box off;    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'ZAxis'), 'FontSize', 5);
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel('Distance', 'FontSize', 6);
    ylabel('Empirical Density', 'FontSize', 6);    
    title('Distribution of Distances to k^0', 'FontSize', 6', 'FontWeight', 'normal');
    
    left_margin = 0.7;
    right_margin = 0.1;    
    top_margin = 0.4;
    bottom_margin = 0.6;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
    
    lgd = legend('FontSize', 5, 'Location', 'NorthEast');
    lgd.Units = 'centimeters';
    lgd.Position(1) = lgd.Position(1)+0.3;
    
    legend('boxoff');
end
