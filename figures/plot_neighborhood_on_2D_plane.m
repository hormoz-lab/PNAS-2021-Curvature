function plot_neighborhood_on_2D_plane(x, y, data, candid_points, candid_r, target_points, x_label, y_label, title_str, series_names)
   
    leg_height = 0.8;

    fig_width = 3.6;
    fig_height = 3.7+leg_height;

    N_series = size(candid_r,2);
    assert(length(data)==N_series);
    
    idx = knnsearch([x(candid_points) y(candid_points)], target_points);
    candid_points = candid_points(idx);
    candid_r = candid_r(idx,:);
    
    for i = 1:N_series
        pw_dist = pdist2(data{i}(candid_points,:), data{i});
        [~, c] = find(pw_dist < candid_r(:,i));
        v{i} = unique(c);
    end
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);
    cols = ['m' 'b' 'g' 'r' 'k'];
    
    h = scatter(x, y, 1, [0.8 0.8 0.8], 'filled', 'DisplayName', '');
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    hold on;
    for i = 1:N_series        
        pl{i} = scatter(x(v{i}), y(v{i}), 1, cols(i), 'filled', 'DisplayName', series_names{i});        
    end
    pl{N_series+1} = scatter(x(candid_points), y(candid_points), 3, 'k', 'filled', 'DisplayName', 'Centre');
    hold off;
    xlim([0 pi]);
    ylim([0 2*pi]);    
    box on;
    
    xticks([0 pi/2 pi]);    
    yticks([0 pi/2 pi 3*pi/2 2*pi]);
    set(gca, 'XTickLabel', {'0'; '$\frac{\pi}{2}$'; '$\pi$'}, ...
             'YTickLabel', {'0'; '$\frac{\pi}{2}$'; '$\pi$'; '$\frac{3\pi}{2}$'; '2$\pi$'}, 'TickLabelInterpreter', 'latex');
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);

    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel(x_label, 'FontSize', 6, 'interpreter', 'latex');
    ylabel(y_label, 'FontSize', 6, 'interpreter', 'latex');    
    t = title(title_str, 'FontSize', 6', 'FontWeight', 'normal');
        
    left_margin = 0.6;
    right_margin = 0.2;    
    top_margin = 0.7;    
    bottom_margin = 0.6;
            
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin+leg_height, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin-leg_height]);
    
    lgd = legend(flipud(vertcat(pl{:})), 'FontSize', 5, 'Location', 'SouthwestOutside');
    if (N_series > 1)
        lgd.NumColumns = 2;
        lgd.Units = 'centimeters';
        lgd.Position = [-0.25, 0, fig_width, leg_height];
    else
        lgd.Units = 'centimeters';
        lgd.Position = [-1.0, 0, fig_width, leg_height];
    end

    legend boxoff;
    
end