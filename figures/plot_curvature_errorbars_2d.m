function sp = plot_curvature_errorbars_2d(S, dS, ref_val, min_pct, max_pct, dS_pct, mask, title_str, which)

    if (nargin == 6)
        full_figure = true;
    else
        full_figure = false;
    end
    
    if (full_figure)
        fig_width = 3.1;
        fig_height = 3.8;        
        figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
               'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
        sp = subplot(1,1,1);
    else
        sp = subplot(1, 4, which);
    end
      
    S = S - ref_val;
    min_val = prctile(S, min_pct);
    max_val = prctile(S, max_pct);
    max_dS_val = prctile(dS, dS_pct);

    ub = linspace(0, max(max_val, max_dS_val), 100);
    lb = linspace(0, max(-min_val, max_dS_val), 100);

    hold on;       
    if (full_figure)
        h = scatter(S+ref_val, dS, 1, 'b', 'filled');
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');    
    else    
        h = scatter(S(mask)+ref_val, dS(mask), 1, 'b', 'filled');
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        h = scatter(S(~mask)+ref_val, dS(~mask), 1, 'g', 'filled');
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    
    h = line( ub+ref_val, ub/2, 'Color', 'r', 'LineWidth', 1.0);
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    line(-lb+ref_val, lb/2, 'Color', 'r', 'DisplayName', sprintf('%d in 95%% CI', ref_val), 'LineWidth', 1.0);
    hold off;
       
    xlim([-lb(end) ub(end)]+ref_val);    
    ylim([0 min(max_dS_val*1.25, max([lb(end)/2; ub(end)/2]))]);
    box off;   
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
        
    if (full_figure)
        xlabel('Scalar Curvature, S', 'FontSize', 6);
    else
        xlabel('S', 'FontSize', 6);
    end
    if (full_figure || which == 1)
        ylabel('Standard Error, \sigma_S', 'FontSize', 6);
    end
    
    if (full_figure)
        lgd = legend('FontSize', 5, 'Location', 'North');
        legend boxoff;
        top_margin = 0.4;
        left_margin = 0.8;
        right_margin = 0.1;        
        bottom_margin = 0.6;
        lgd_width = 2.5;
        lgd_height = 0.3;        
        set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
        lgd.Units = 'centimeters';
        lgd.Position = [fig_width-lgd_width fig_height-lgd_height lgd_width lgd_height];        
    else
        title({title_str, ''}, 'FontSize', 6, 'FontWeight', 'normal', 'interpreter', 'latex');        
    end
end
