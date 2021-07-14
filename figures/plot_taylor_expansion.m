function plot_taylor_expansion(x1, x2, score, analytic, center_eps)

    fig_width = 5.4;
    fig_height = 3.8;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);   
    hold on;
    colormap(jet)
    surf(x1, x2, analytic, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    scatter3(score(:,1),score(:,2),score(:,3),30,score(:,3),'filled');    
    view([45 30]);    
    hold off;
    
    xlim([-center_eps center_eps]);
    ylim([-center_eps center_eps]); 
    zticks([-0.01 0.01]);
    zticklabels({'-0.1'; '0'});
    
    box off;    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'ZAxis'), 'FontSize', 5);
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel('Tangent $t_1$', 'FontSize', 6, 'interpreter', 'latex');
    ylabel('Tangent $t_2$', 'FontSize', 6, 'interpreter', 'latex');
    zlabel('Normal $n_1$', 'FontSize', 6, 'interpreter', 'latex');    
    title('$n_k \approx C_k+\sum_{i,j} h_{ij}^k t_i t_j$', 'FontSize', 8, 'interpreter', 'latex');
        
    left_margin = 0.9;
    right_margin = 0.1;    
    top_margin = 0.5;
    bottom_margin = 0.5;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
   
end