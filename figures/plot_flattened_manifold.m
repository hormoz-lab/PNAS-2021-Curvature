function plot_flattened_manifold(X, C, title_str, xlab, ylab, xl, yl, xt, yt, xtl, ytl)

    if (nargin < 10)
        xtl = cellstr(num2str(xt'));
        ytl = cellstr(num2str(yt'));
    end

    fig_width = 2.5;
    fig_height = 3.0;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);
    scatter(X(:,1), X(:,2), 1, C,'filled');        
    colormap(jet);
    xticks(xt);
    yticks(yt);    
    set(gca, 'XTickLabel', xtl, ...
             'YTickLabel', ytl, 'TickLabelInterpreter', 'latex');
    xlim(xl);
    ylim(yl);
    
    box off;
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    
    xlabel(xlab, 'FontSize', 6);
    ylabel(ylab, 'FontSize', 6);    
    title(title_str, 'FontSize', 6', 'FontWeight', 'normal');
    
    left_margin = 0.8;
    right_margin = 0.1;    
    top_margin = 0.3;
    bottom_margin = 0.6;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
    
end
   