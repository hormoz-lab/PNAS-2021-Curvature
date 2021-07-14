function plot_trace_fitted_curvature(x, S, which_trace)

    assert(length(x) == size(S,1));
    assert(length(x) == size(S,2));
    
    fig_width = 2.6;
    fig_height = 3.2;
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);
    hold on;
    h = pcolor(x, x, S');    
    set(h, 'EdgeColor', 'none');
    colormap(jet);    
    
    cb_lo = nanmin(S(:));
    cb_hi = nanmax(S(:));
    if (cb_hi>10)
        round_fact = 0;
    else
        round_fact = 1;
    end
    if (isnan(cb_lo))
        cb_lo = 1.5;
    end
    if (isnan(cb_hi))
        cb_hi = 2.5;
    end
        
    cb = colorbar('FontSize', 5, 'Location', 'SouthOutside');    
    caxis([cb_lo cb_hi]);
    cb.Ticks = [cb_lo cb_hi];
    cb.TickLabels = {num2str(round(cb_lo,round_fact)); num2str(round(cb_hi,round_fact))};    
    cb.AxisLocation = 'out';    
    patch([0 x(end) x(end)], [0 0 x(end)], 'black', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    axis tight;    
    box on;
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel('x_1', 'FontSize', 6);
    ylabel('x_2', 'FontSize', 6);    
    title(sprintf('S from %s', which_trace), 'FontSize', 6', 'FontWeight', 'normal', 'interpreter', 'latex');
    
    left_margin = 0.6;
    right_margin = 0.0;    
    top_margin = 0.3;
    bottom_margin = 1.3;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
    
    cb.Units = 'centimeters';
    cb.Position = [left_margin 3*bottom_margin/8 fig_width-left_margin-right_margin bottom_margin/8];
    yl = ylabel(cb, 'Scalar Curvature, S');
    p = get(yl,'position');
    p(2) = 0.3*p(2);
    set(yl,'position',p);
    
end