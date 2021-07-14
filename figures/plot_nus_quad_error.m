function plot_nus_quad_error(dat, title_str, N1, N2, min_pct, max_pct, dS_pct)

    fig_width = 8.8;
    fig_height = 3.8;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
       
    sp_title = {'$\sigma_{\overline{h}}$'; '$r_{5}$'; '$r_{50}$'; '$r_{95}$'};
    
    for i = 1:4    
        sp{i} = plot_curvature_errorbars_2d(dat{i}.S, dat{i}.dS, 0, min_pct(i), max_pct(i), dS_pct(i), [true(N1,1); false(N2,1)], sp_title{i}, i);
    end
    
    top_margin = 1.0;
    bottom_margin = 0.6;
    left_margin = 0.9;
    right_margin = 0.1; 
    tight_margin = 0.5;
    
    panel_width = (fig_width-left_margin-right_margin-3*tight_margin)/4;
    panel_height = fig_height-bottom_margin-top_margin;
    
    set(sp{1}, 'Units', 'centimeters', 'Position', [left_margin + 0*(tight_margin+panel_width), bottom_margin, panel_width, panel_height]);
    set(sp{2}, 'Units', 'centimeters', 'Position', [left_margin + 1*(tight_margin+panel_width), bottom_margin, panel_width, panel_height]);
    set(sp{3}, 'Units', 'centimeters', 'Position', [left_margin + 2*(tight_margin+panel_width), bottom_margin, panel_width, panel_height]);
    set(sp{4}, 'Units', 'centimeters', 'Position', [left_margin + 3*(tight_margin+panel_width), bottom_margin, panel_width, panel_height]);
    
    h = suptitle(title_str);
    h.FontSize = 6;
    h.FontWeight = 'normal';
    h.Interpreter = 'tex';
    
end