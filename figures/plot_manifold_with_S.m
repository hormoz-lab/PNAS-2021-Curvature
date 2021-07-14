function plot_manifold_with_S(X, S, clims, title_str, cb_label)

    if (nargin == 4)
        cb_label = 'Scalar Curvature, S';
    end
        
    fig_width = 3.1;
    fig_height = 3.8;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);
    scatter3(X(:,1),X(:,2), X(:,3), 1, S,'filled');
    view([45 30])
    colormap(jet);
    caxis(round(clims,1));        
    cb = colorbar('FontSize', 5, 'Location', 'SouthOutside', 'TickDirection', 'out');    
    cb.Ticks = round(clims,1);
    axis tight;
    box off;
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'ZAxis'), 'FontSize', 5);
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel('x', 'FontSize', 6);
    ylabel('y', 'FontSize', 6);
    zlabel('z', 'FontSize', 6);
    title(title_str, 'FontSize', 6', 'FontWeight', 'normal', 'interpreter', 'latex');
    
    left_margin = 0.7;
    right_margin = 0.2;    
    top_margin = 0.4;
    bottom_margin = 1.0;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
   
    cb.Units = 'centimeters';
    cb.Position = [left_margin 3*bottom_margin/8 fig_width-left_margin-right_margin bottom_margin/8];
    yl = ylabel(cb, cb_label);
    p = get(yl,'position');
    p(2) = 0.3*p(2);
    set(yl,'position',p);
    
end