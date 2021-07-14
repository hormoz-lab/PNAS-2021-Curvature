function plot_SS_stat(pos, neg, which_SS, x_label, y_label, title_str)

    fig_width = 2.5;
    fig_height = 3.8;
    effect_size = 50;
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);
    hold on;
    b = bar([1:length(which_SS)], [pos' mean(neg,2)], 'EdgeColor', 'none');
    errorbar([1:length(which_SS)]+0.15, mean(neg,2), std(neg,[], 2), 'Marker', 'none', 'LineWidth', 1.0', 'Color', 'black', 'LineStyle', 'none');
    sig_K = find(pos'>max(neg, [], 2) & pos'>effect_size);
    plot(sig_K, pos(sig_K)+0.1*max(pos(sig_K)), 'Color', 'black', 'Marker', '*', 'LineStyle', 'none', 'MarkerSize', 5);
    hold off;
        
    axis tight;
    yl = ylim;
    ylim([yl(1) yl(2)*1.05]);
    
    box off;    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xticks([1:length(which_SS)+1]);
    xticklabels([cellstr(num2str(which_SS')); 'Rand']);
    xlabel(x_label, 'FontSize', 6);
    ylabel(y_label, 'FontSize', 6);    
    title(title_str, 'FontSize', 6, 'FontWeight', 'normal');
    
    left_margin = 0.7;
    right_margin = 0.1;    
    top_margin = 0.7;
    bottom_margin = 1.2;
            
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
    
    leg_height = 0.4;
    lgd = legend({'Matched'; 'Random'}, 'Location', 'SouthOutside', 'FontSize', 5, 'NumColumns', 1);
    legend('boxoff');
    lgd.Units = 'Centimeters';
    leg_pos = lgd.Position;
    lgd.Position = [leg_pos(1)+0.2 0.1 leg_pos(3) leg_height];
    

end