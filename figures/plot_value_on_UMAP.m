function plot_value_on_UMAP(x, y, S, title_str, clims)

    fig_width = 4.3;
    fig_height = 3.8;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    sp1 = subplot(1,1,1);
    colormap jet; scatter(x, y, 1, S, 'filled');     
    if (isscalar(clims))
        if (clims==true)
            caxis([prctile(S,1) prctile(S,99)]);
        else
            caxis([nanmin(S) nanmax(S)]);
        end
    else    
        caxis(round(clims,1));        
    end    
    cb = colorbar('FontSize', 5, 'Location', 'EastOutside');
    
    axis tight;
    axis off;
    box off;
    t = title(title_str, 'FontSize', 6', 'FontWeight', 'normal');
        
    left_margin = 0.30;
    right_margin =1.2;    
    top_margin = 0.4;
    bottom_margin = 0.1;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);
    
    if (iscell(title_str))
       t.Position(2) = t.Position(2)-1.5; 
    end

end