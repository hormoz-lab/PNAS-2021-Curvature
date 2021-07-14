function plot_stat_histogram(X, x_label, do_clip, vbar)

    fig_width = 4.25;
    fig_height = 3.8;
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);
    if (do_clip)
        X = X(X>prctile(X,1) & X<prctile(X,99));
    end
    [freqs, xvals] = histcounts(X,100);
    freqs = 100*freqs/sum(freqs);
    xvals = movmean(xvals, 2, 'Endpoints', 'discard');
    bar(xvals, freqs, 'EdgeColor', 'none');
    
    if (nargin == 4)
        hold on;
        for i = 1:length(vbar)
            line([vbar(i) vbar(i)], [0, max(freqs)], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.0);
        end
        hold off;
    end
    axis tight;
    
    box off;    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel(x_label, 'FontSize', 6);
    ylabel('% Points', 'FontSize', 6);
    
    left_margin = 0.7;
    right_margin = 0.1;    
    top_margin = 0.1;
    bottom_margin = 0.6;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);

end