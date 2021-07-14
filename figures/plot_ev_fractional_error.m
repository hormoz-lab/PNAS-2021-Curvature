function plot_ev_fractional_error(frac_error, series_names, f)
    
    assert(length(series_names)==size(frac_error,2));
    n = size(frac_error,1);

    fig_width = 8.8;
    fig_height = 4;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,2,1);
    hold on;
    plot(1:n, frac_error(:,1), 'LineWidth', 1.0, 'Color', 'red', 'DisplayName', series_names{1});
    plot(1:n, frac_error(:,2), 'LineWidth', 1.0, 'Color', 'blue', 'DisplayName',  series_names{2});
    plot(1:n, frac_error(:,3), 'LineWidth', 1.0, 'Color', 'cyan', 'DisplayName',  series_names{3});
    plot(1:n, frac_error(:,4), 'LineWidth', 1.0, 'Color', 'black', 'DisplayName',  series_names{4});    
    hold off;    
    axis tight;
    xlim([1, floor(sqrt(n))].^2);    
    xtickvals = [1:floor(sqrt(n))].^2;    
    xticks(xtickvals);
    xticklabels(cellstr(num2str(xtickvals')));
    yl = ylim;
    ylim([yl(1) 1]);
    box off;
        
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel('Eigenvalue', 'FontSize', 6);
    ylabel('Fractional Error', 'FontSize', 6);
    title('Estimator Eigenspectra', 'FontSize', 6', 'FontWeight', 'normal');
    
    lgd1 = legend('FontSize', 5, 'Location', 'Southeast', 'NumColumns', 2);
    legend('boxoff');
    
    sp2 = subplot(1,2,2);    
    ax = gca;
    set(ax, 'yaxislocation', 'right');
    hold on;    
    n_grid = linspace(4, 10, 100)';
    N = 10^4;
    prefactor = mean(frac_error(37:49,1)*N^0.25/(log(N))^(3/8));
    error_curve = arrayfun(@(n) prefactor*log(10^n)^(3/8)/(10^n)^(1/4), n_grid);
    plot(n_grid, error_curve, 'LineWidth', 1.0, 'Color', 'red');
    line(xlim,[frac_error(end,1)*f frac_error(end,1)*f],'Color','green', 'LineWidth', 1.0, 'LineStyle', '--')
    hold off;
        
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel('log_{10} Sample Size', 'FontSize', 6);
    ylabel('Fractional Error', 'FontSize', 6);    
    title('Best Case \lambda_{37},...,\lambda_{49} Error', 'FontSize', 6', 'FontWeight', 'normal');
     
    lgd2 = legend({series_names{1}; 'Target Error'}, 'Location', 'Northeast', 'FontSize', 5);
    legend('boxoff');
    
    left_margin = 0.8;
    right_margin = 0.8;    
    mid_margin = 0.4;
    top_margin = 0.4;
    bottom_margin = 1.25;
    
    panel_width = (fig_width-left_margin-right_margin-mid_margin)/2;
    
    linkaxes([sp1 sp2], 'y');
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, panel_width, fig_height-bottom_margin-top_margin]);
    set(sp2, 'Units', 'centimeters', 'Position', [left_margin+panel_width+mid_margin, bottom_margin, panel_width, fig_height-bottom_margin-top_margin]);
    
    lgd1.Units = 'centimeters';
    lgd_height = 0.6;
    lgd1.Position = [0 0 fig_width lgd_height];
             
end