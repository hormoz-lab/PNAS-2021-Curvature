function plot_Sn_distribution(n, S)

    assert(length(n)==length(S));

    fig_width = 6;
    fig_height = 3.6;
        
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
       
    sp1 = subplot(1,1,1);
    hold on;
    max_f = 0;
    min_S = min(cellfun(@(x) nanmin(x(x>0)), S));
    max_S = max(cellfun(@(x) nanmax(x(x>0)), S));
    
    for i = 1:length(n)       
        histogram(log2(S{i}(S{i}>0)),linspace(log2(min_S), log2(max_S), 201),'Normalization','pdf');    
        [f,x] = ksdensity(log2(S{i}(S{i}>0)));
        plot(x,f,'LineWidth',1,'Color','k');
        max_f = max(max(f),max_f);
    end
    for i = 1:length(n)        
        line(log2(n(i)*(n(i)-1))*[1 1], [0 max_f],'Color','r','LineStyle', '--', 'LineWidth', 1.0);
        text(log2(n(i)*(n(i)-1))-0.5, max_f-1, sprintf('d=%d', n(i)), 'FontSize', 5);
    end
    hold off;
   
    xlim([log2(1) log2(64)]);
    set(gca,'xtick', 0:8)
    set(gca,'xticklabels',2.^[0:8])
    ylim([0 max_f]);
    yticks([]);
       
    box off;    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'ZAxis'), 'FontSize', 5);
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel('Scalar Curvature, S', 'FontSize', 6);
    ylabel('Empirical Density', 'FontSize', 6);
    title('10K Points from $\mathcal{S}^d$', 'FontSize', 6', 'FontWeight', 'normal', 'interpreter', 'latex');
    
    left_margin = 0.5;
    right_margin = 0.1;    
    top_margin = 0.3;
    bottom_margin = 0.6;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);

end
