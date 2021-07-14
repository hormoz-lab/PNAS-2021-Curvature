function plot_downsample_correlation(x, y, x_label, y_label, title_str, corr_type, mega, clip)

    assert(length(x) == length(y));
    
    mask = ~(isnan(x) | isnan(y));
    x = x(mask);
    y = y(mask);
    
    n = floor(log10(length(x)));
    
    [R, p] = corr(x, y, 'Type', corr_type, 'Rows', 'complete');
    if (isequal(corr_type, 'Spearman'))
        [~, x] = ismember(x, sort(x));
        [~, y] = ismember(y, sort(y));
    end

    if (mega)
        fig_width = 14.3;
        fig_height = 13.8;
        fs_sm = 12;
        fs_lg = 14;   
        mk_sz = 8;
    else
        fig_width = 3.75;
        fig_height = 3.8;
        fs_sm = 5;
        fs_lg = 6;
        mk_sz = 2;
    end

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);
    scatter(x, y, mk_sz/n^2, 'b', 'filled');     
        
    axis tight;
    if (nargin == 8)
        xlim([prctile(x, clip(1)), prctile(x, clip(2))]);
        ylim([prctile(x, clip(1)), prctile(x, clip(2))]);
    end
    
    box off;    
    set(get(gca, 'XAxis'), 'FontSize', fs_sm);
    set(get(gca, 'YAxis'), 'FontSize', fs_sm);    
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');    
    xlabel(x_label, 'FontSize', fs_lg);
    ylabel(y_label, 'FontSize', fs_lg);
    
    if (isequal(corr_type, 'Spearman'))
        ax = gca;    
        ax.XAxis.Exponent = n;
        ax.YAxis.Exponent = n;
    end

    title(title_str, 'FontSize', fs_lg, 'FontWeight', 'normal');
            
    if (mega)
        left_margin = 1.7;
        right_margin = 0.1;
        top_margin = 1.3;
        bottom_margin = 1.7;
    else
        left_margin = 0.7;
        right_margin = 0.2;
        top_margin = 0.7;
        bottom_margin = 0.7;
    end
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, ...
             fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
    if (p<1e-6)
        str = {sprintf('\\rho=%.2f', R); 'p<1.00e-6'};
    else
        str = {sprintf('\\rho=%.2f', R); sprintf('p=%.2e', p)};
    end
         
    if (mega)
       box_pos = [0.75 0.06 0.2 0.2];
    else
        box_pos = [0.6 0.2 0.2 0.2];
    end
    
    an = annotation('textbox', box_pos, 'String',str,'FitBoxToText','on');
    an.FontSize = fs_sm;
    an.FontWeight = 'normal';
    an.BackgroundColor = [1 1 1];
    

end