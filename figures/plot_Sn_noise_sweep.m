function plot_Sn_noise_sweep(dat, noise_levels, dims)
    
    assert(size(dat,1) == length(dims));
    assert(size(dat,2) == length(noise_levels));
    assert(length(noise_levels) == 3);
    
    fig_width = 9.0;
    fig_height = 4;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
       
    for j = 1:3
        sp{j} = subplot(1,3,j);
        hold on;
        for i = 1:length(dims)
            [f,x] = ksdensity(dat{i,j}.S);
            plot(x, f, 'LineWidth',1);
        end        
        axis tight;
        yl = ylim;
        h = line([2 2], yl, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.0);
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
        hold off;    
        box off;    
        set(gca, 'ytick', []);
        set(get(gca, 'XAxis'), 'FontSize', 5);
        set(get(gca, 'YAxis'), 'FontSize', 5);
        set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
        if (j == 2)
            xlabel('Scalar Curvature, S', 'FontSize', 6);
        end
        if (j == 1)
            ylabel('Empirical Density', 'FontSize', 6);
        end
        title(sprintf('\\sigma=%.2f', noise_levels(j)), 'FontSize', 6, 'FontWeight', 'normal');
    end
    
    lgd = legend(cellstr(num2str(dims')), 'FontSize', 5, 'Location', 'SouthOutside', 'NumColumns', 4);
    title(lgd, sprintf('Ambient Dimension, n'), 'FontSize', 6, 'FontWeight', 'normal');
    legend('boxoff');
    
    left_margin = 0.5;
    right_margin = 0.1;
    mid_margin = 0.1;
    top_margin = 0.4;
    bottom_margin = 1.6;
    
    panel_width = (fig_width-left_margin-right_margin-2*mid_margin)/3;
    
    set(sp{1}, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, panel_width, fig_height-bottom_margin-top_margin]);
    set(sp{2}, 'Units', 'centimeters', 'Position', [left_margin+mid_margin+panel_width, bottom_margin, panel_width, fig_height-bottom_margin-top_margin]);
    set(sp{3}, 'Units', 'centimeters', 'Position', [left_margin+2*mid_margin+2*panel_width, bottom_margin, panel_width, fig_height-bottom_margin-top_margin]);
    
    lgd.Units = 'centimeters';
    lgd_height = 1.0;
    lgd.Position = [0 0 fig_width lgd_height];    
    
end
