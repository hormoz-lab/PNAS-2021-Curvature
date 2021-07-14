function plot_neighborhood_on_UMAP(x, y, data, candid_points, candid_r, target_points)
    
    idx = knnsearch([x(candid_points) y(candid_points)], target_points);
    candid_points = candid_points(idx);
    candid_r = candid_r(idx);
    
    pw_dist = pdist2(data(candid_points,:), data);
    [~, c] = find(pw_dist < candid_r);
    v = unique(c);
    
    fig_width = 3.35;
    fig_height = 3.8;
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);
    
    h = scatter(x, y, 1, [0.8 0.8 0.8], 'filled', 'DisplayName', '');
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    hold on;    
    pl{1} = scatter(x(v), y(v), 1, 'm', 'filled', 'DisplayName', 'Neighborhood');    
    pl{2} = scatter(x(candid_points), y(candid_points), 3, 'k', 'filled', 'DisplayName', 'Centre');
    hold off;

    axis tight;    
    box off;
    axis off;
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);

    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    
    lgd = legend(flipud(vertcat(pl{:})), 'FontSize',5, 'Location', 'North', 'NumColumns', 1);
    
    leg_width = 2;
    leg_height = 0.5;
    
    left_margin = 0.1;
    right_margin = 0.1;    
    top_margin = 0.4;
    bottom_margin = 0.1;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
    
    lgd.Units = 'centimeters';
    lgd.Position = [0, fig_height-leg_height, leg_width, leg_height];
    legend boxoff;
    
end