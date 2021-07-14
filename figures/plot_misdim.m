function plot_misdim(dat, true_S, true_dim, title_str)

    fig_width = 8.7;
    fig_height = 4;
        
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);  
    
    c = ['c', 'm', 'g'];
        
    transform_log = @(x) sign(x).*log2(1+abs(x));
    transform_nat = @(x) sign(x).*(2.^abs(x)-1);
    
    sp1 = subplot(1,3,1);
    mask = cellfun(@(x) x.S>=prctile(x.S,5) & x.S<=prctile(x.S,95), dat, 'un', false);
    S_min = min(cellfun(@(x,m) nanmin(x.S(m)), dat, mask));
    S_max = max(cellfun(@(x,m) nanmax(x.S(m)), dat, mask));
    S_min = transform_log(S_min);
    S_max = transform_log(S_max);
    hold on;
    for i = 1:length(dat)
        temp = transform_log(dat{i}.S(mask{i}));     
        xedges = linspace(S_min, S_max, 50);
        counts = histcounts(temp, xedges);
        bar(movmean(xedges,2, 'endpoints', 'discard'), counts, 'FaceColor', c(i), 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'BarWidth', 1.0);
    end
    axis tight;
    yl = ylim;
    line(transform_log(true_S)*[1 1], yl, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.0);
    hold off;        
    box off; 
    
    set(gca, 'ytick', []);
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
      
    xt_range = transform_log([-64 -16 -4 -1 0 1 4 16 64]);    
    xt = xt_range(find(xt_range<=S_min, 1, 'last'):find(xt_range>=S_max, 1, 'first'));
    xticks(xt);
    xticklabels(cellstr(num2str([transform_nat(xt)]')));    
    
    xlabel('Scalar Curvature, S', 'FontSize', 6);    
    ylabel('Empirical Density', 'FontSize', 6);
    
    sp2 = subplot(1,3,2);
    mask = cellfun(@(x) x.dS>=prctile(x.dS,5) & x.dS<=prctile(x.dS,95), dat, 'un', false);
    dS_min = min(cellfun(@(x,m) nanmin(x.dS(m)), dat, mask));
    dS_max = max(cellfun(@(x,m) nanmax(x.dS(m)), dat, mask));
    dS_min = log2(dS_min);
    dS_max = log2(dS_max);    
    hold on;
    for i = 1:length(dat)
        xedges = linspace(dS_min, dS_max, 50);
        counts = histcounts(log2(dat{i}.dS(mask{i})), xedges);
        bar(movmean(xedges,2, 'endpoints', 'discard'), counts, 'FaceColor', c(i), 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'BarWidth', 1.0);
    end
    axis tight;    
    hold off;
    
    box off;    
    set(gca, 'ytick', []);
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    
    xlabel('Standard Error, log_2 \sigma_S', 'FontSize', 6);
        
    sp3 = subplot(1,3,3);    
    hold on;
    for i = 1:length(dat)
        xedges = linspace(0, 1, 10);
        counts = histcounts(dat{i}.gof, xedges);
        plot(movmean(xedges,2, 'endpoints', 'discard'), counts, 'LineWidth', 1.0, 'Color', c(i));        
    end
    axis tight;    
    hold off;
    
    box off;
    set(gca, 'ytick', []);
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel('GOF p-value', 'FontSize', 6);
        
    h = suptitle(title_str);
    set(h, 'FontSize', 6, 'FontWeight', 'normal', 'interpreter', 'latex');
    
    lgd_entries = cellfun(@(x,y) [num2str(mean(x.ball_r), '%#-.2f') ' (' num2str(y, '%1d') ')'], ...
                                  dat, num2cell([true_dim-1, true_dim+1, true_dim]), 'un', false);
    
    lgd = legend(lgd_entries, 'FontSize', 5, 'Location', 'SouthOutside', 'NumColumns', 3);
    title(lgd, 'Average Ball Radius, r (Manifold Dimension, d)', 'FontSize', 6, 'FontWeight', 'normal');
    legend('boxoff');
        
    left_margin = 0.4;
    right_margin = 0.1;
    mid_margin = 0.1;
    top_margin = 0.4;
    bottom_margin = 1.4;
    
    panel_width = (fig_width-left_margin-right_margin-2*mid_margin)/3;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, panel_width, fig_height-bottom_margin-top_margin]);
    set(sp2, 'Units', 'centimeters', 'Position', [left_margin+mid_margin+panel_width, bottom_margin, panel_width, fig_height-bottom_margin-top_margin]);
    set(sp3, 'Units', 'centimeters', 'Position', [left_margin+2*mid_margin+2*panel_width, bottom_margin, panel_width, fig_height-bottom_margin-top_margin]);
    
    lgd.Units = 'centimeters';
    lgd_height = 0.8;
    lgd.Position = [0 0 fig_width lgd_height];
    
end