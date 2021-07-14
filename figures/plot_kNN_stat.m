function plot_kNN_stat(pos, neg, which_K, x_label, y_label, title_str)

    fig_width = 4.25;
    fig_height = 3.8;
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);
    hold on;
    bar([1:length(which_K)], pos, 'EdgeColor', 'none', 'FaceColor', 'b');
    bar(length(which_K)+1, mean(neg), 'EdgeColor', 'none', 'FaceColor', 'r');
    errorbar(length(which_K)+1, mean(neg), std(neg), 'Marker', 'none', 'LineWidth', 1.0', 'Color', 'black');
    sig_K = find(pos>max(neg));
    plot(sig_K, pos(sig_K)+0.1*max(pos(sig_K)), 'Color', 'black', 'Marker', '*', 'LineStyle', 'none', 'MarkerSize', 5);
    hold off;
    
    axis tight;
    yl = ylim;
    ylim([0 yl(2)*1.05]);
    
    box off;    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xticks([1:length(which_K)+1]);
    xticklabels([cellstr(num2str(which_K')); 'Rand']);
    xlabel(x_label, 'FontSize', 6);
    ylabel(y_label, 'FontSize', 6);    
    title(title_str, 'FontSize', 6, 'FontWeight', 'normal');
    
    left_margin = 0.7;
    right_margin = 0.1;    
    top_margin = 0.3;
    bottom_margin = 0.6;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);

end