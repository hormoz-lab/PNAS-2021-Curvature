function plot_S_on_2D_plane(x, y, S, x_label, y_label, title_str, do_clip, tall)

    if (nargin == 7)
        tall = false;
    end
    
    fig_width = 3.5;
    if (tall)
        fig_height = 3.7;
    else        
        fig_height = 3.4;
    end

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);
    colormap jet; scatter(x, y, 1, S, 'filled'); 
    cb = colorbar('FontSize', 5, 'Location', 'EastOutside');

    if (do_clip)
        caxis([prctile(S,5) prctile(S,95)]);
    else
        caxis([nanmin(S) nanmax(S)]);
    end
    xlim([0 pi]);
    ylim([0 2*pi]);    
    box on;
    xticks([0 pi/2 pi]);    
    yticks([0 pi/2 pi 3*pi/2 2*pi]);
    set(gca, 'XTickLabel', {'0'; '$\frac{\pi}{2}$'; '$\pi$'}, ...
             'YTickLabel', {'0'; '$\frac{\pi}{2}$'; '$\pi$'; '$\frac{3\pi}{2}$'; '2$\pi$'}, 'TickLabelInterpreter', 'latex');
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'ZAxis'), 'FontSize', 5);
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel(x_label, 'FontSize', 6, 'interpreter', 'latex');
    ylabel(y_label, 'FontSize', 6, 'interpreter', 'latex');    
    t = title(title_str, 'FontSize', 6', 'FontWeight', 'normal');
    
    left_margin = 0.6;
    right_margin = 0.9;    
    if (tall)
        top_margin = 0.7;
    else
        top_margin = 0.4;
    end
    bottom_margin = 0.6;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
   
    annotation('textbox', 'Position', [cb.Position(1) cb.Position(2)-0.05 0.05 0.05], 'string', 'S', 'EdgeColor', 'none', 'FontSize', 6, 'FontWeight', 'normal');
    
end