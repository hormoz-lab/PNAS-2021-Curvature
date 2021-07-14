function plot_heat_traces(x, traces, series_name)

    assert(length(traces) == 6);
    assert(length(series_name)==6);
    assert(all(length(x)==cellfun(@length, traces)));

    fig_width = 8.8;
    fig_height = 4;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);
    hold on;
    plot(x, traces{1}, 'LineWidth', 2.5, 'Color', 'black', 'DisplayName', series_name{1});
    plot(x, traces{2}, 'LineWidth', 1.0, 'Color', 'blue',  'DisplayName', series_name{2});
    plot(x, traces{3}, 'LineWidth', 1.0, 'Color', 'red',   'DisplayName', series_name{3});
    plot(x, traces{4}, 'LineWidth', 1.0, 'Color', 'blue',  'DisplayName', series_name{4}, 'LineStyle', '--');
    plot(x, traces{5}, 'LineWidth', 1.0, 'Color', 'red',   'DisplayName', series_name{5}, 'LineStyle', '--');
    plot(x, traces{6}, 'LineWidth', 1.0, 'Color', 'green', 'DisplayName', series_name{6});
    hold off;
    x_lb = 0.3;
    x_ub = 1.05;
    max_z = 30;
    patch([0 0 x_lb x_lb], [0 max_z max_z 0], 'black', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    patch([x_ub x_ub x(end) x(end)], [0 max_z max_z 0], 'black', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    hold off;
    axis tight;
    ylim([0 max_z]);
    xlim([x(1) x(end)]);
    xtickvals = [x(1) x_lb x_ub x(end)];
    xticks(xtickvals);
    xticklabels({num2str(x(1)); 'x_1'; 'x_2'; num2str(x(end))});
    box off;
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'ZAxis'), 'FontSize', 5);
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel('x', 'FontSize', 6);
    ylabel('Heat-Trace', 'FontSize', 6);
   
    left_margin = 0.7;
    right_margin = 0.2;    
    top_margin = 0.1;
    bottom_margin = 0.75;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
    
    lgd = legend(series_name, 'FontSize', 5, 'interpreter', 'Latex');
    legend('boxoff');
   
    lgd_width = 4.4;
    lgd_height = 1;
    lgd.NumColumns = 2;
    lgd.Units = 'centimeters';    
    lgd.Position = [fig_width-right_margin-lgd_width, bottom_margin, lgd_width, lgd_height];
            
end