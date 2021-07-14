function plot_curvature_with_errorbars(analytic, S, dS, x_label, x_coords)

    fig_width = 3.1;
    fig_height = 3.6;   
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
  
    mask = ~isnan(S);    
    analytic = analytic(mask);    
    S = S(mask);
    dS = dS(mask);    
    x_coords = x_coords(mask);
            
    sp1 = subplot(1,1,1);
    hold on;    
    plot(x_coords, S+2*dS, 'LineWidth', 0.5, 'color', [0.5 0.5 0.5]);
    plot(x_coords, S-2*dS, 'LineWidth', 0.5, 'color', [0.5 0.5 0.5]);
    plot(x_coords, S, 'o', 'MarkerSize', 0.1, 'color', 'black');
    plot(x_coords, analytic, 'LineWidth', 1.0, 'LineStyle', '--', 'color', 'red');
    axis tight;
    ylim([prctile(S-2*dS, 0.5), prctile(S+2*dS, 99.5)]);
    
    box off;    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');    
    xlabel(x_label, 'FontSize', 6);    
    ylabel('Scalar Curvature 95% CI (S\pm2\sigma_S)', 'FontSize', 6);
    
    left_margin = 0.8;
    right_margin = 0.1;    
    top_margin = 0.4;
    bottom_margin = 0.6;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);

end
