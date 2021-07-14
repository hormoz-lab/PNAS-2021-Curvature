function plot_cell_annotation(coords, cell_type)

    unannotated = cellfun(@isempty, cell_type);
    coords(unannotated,:) = [];
    cell_type(unannotated,:) = [];

    [cell_type, ~, which] = unique(cell_type, 'stable');
    [~, cell_order] = sort(accumarray(which, 1), 'descend');   
    N_types = max(which);
    
    fig_width = 3.0;
    fig_height = 3.5;
    if (length(coords)>1e4)
        fig_width = 4.2;
        fig_height = 4.9;
    end
    
    N_cols = ceil(N_types/13);
        
    lgd_width = 4.0*N_cols;
    fig_width = fig_width+lgd_width;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
       
    c1 = [1 0 0;
          0 0 1;
          0 1 0;
          1 1 0;
          0 1 1;
          1 0 1];
    c2 = unique(lines, 'rows', 'stable');
    c3 = [128 128 128;
          233 185 210;
          240 185 140;
          163 181 141;
           83 177 247;
          169 145 145;
          212  94  94;
          152 174 246;
           91 165 100;
          192 192 192]/255;
    
    colors = [c1; c2; c3; 2*c1/3; [1 1 1]/3; 1-(1-c1)/3; [0 0 0]];
            
    sp1 = subplot(1,1,1);
    hold on;
    for i = 1:N_types
        pl{i} = scatter(coords(which==cell_order(i), 1), coords(which==cell_order(i), 2), 1, colors(N_types-i+1,:), 'filled', 'DisplayName', cell_type{cell_order(i)});
    end
    axis tight;
    axis off;
    box off;    
    
    lgd = legend('FontSize', 5, 'Location', 'EastOutside', 'NumColumns', N_cols);
    
    legend('boxoff');
        
    top_margin = 0.1;
    left_margin = 0.0;
    right_margin = 0.0;        
    bottom_margin = 0.1;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin-lgd_width, fig_height-bottom_margin-top_margin]);
    
    lgd.Units = 'centimeters';    
    lgd.Position = [fig_width-lgd_width, bottom_margin, lgd_width, fig_height-bottom_margin];

end