function plot_S2_noise_sweep(S, noise_level)
   
    assert(length(S)==length(noise_level));
    
    fig_width = 4.5;
    fig_height = 4;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,2,1);
    hold on;
    for i = 1:length(noise_level)
        [f,x] = ksdensity(S{i}(~isnan(S{i})));
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
    xlabel('Scalar Curvature, S', 'FontSize', 6);
    ylabel('Empirical Density', 'FontSize', 6);
    title('10K Points from Noisy $\mathcal{S}^2$', 'FontSize', 6, 'FontWeight', 'normal', 'interpreter', 'latex');
    
    lgd = legend(cellstr(num2str(noise_level, '%#-.1g')), 'FontSize', 5, 'Location', 'East');
    title(lgd, '\sigma', 'FontSize', 6, 'FontWeight', 'normal');
    legend('boxoff');
    
    left_margin = 0.5;
    right_margin = 0.1;
    top_margin = 0.4;
    bottom_margin = 0.6;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
    
    lgd.Units = 'centimeters';    
    lgd_width = 1.5;
    lgd.Position = [fig_width-lgd_width-right_margin 2*bottom_margin lgd_width fig_height-top_margin-2*bottom_margin];
    
end
