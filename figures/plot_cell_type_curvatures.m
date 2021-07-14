function plot_cell_type_curvatures(cell_type, S, N_neighbors, method)

    unannotated = cellfun(@isempty, cell_type);
    
    cell_type(unannotated,:) = [];
    S(unannotated,:) = [];
    N_neighbors(unannotated,:) = [];

    [cell_type, ~, which_type] = unique(cell_type, 'stable');    
    [~, cell_order] = sort(accumarray(which_type, 1), 'descend');
    cell_type = cell_type(cell_order);
    [~, which_type] = ismember(which_type, cell_order);
    N_type = length(cell_type);
    
    [r, c] = find(tril(ones(N_type),-1));
            
    if (isequal(method, 'Global'))
        [mu, se, df] = arrayfun(@(i) get_t_prop(S(which_type==i), N_neighbors, method), [1:N_type]');    
    else
        [mu, se, df] = arrayfun(@(i) get_t_prop(S(which_type==i), N_neighbors(which_type==i), method), [1:N_type]');    
    end
    [p, mu_diff, t_stat, df_joint] = arrayfun(@(i,j) welch_t_test(mu(i), mu(j), se(i), se(j), df(i), df(j)), r, c);
   
    p_crit = benjamini_hochberg(p(~isnan(p)), 0.05);
    
    left_margin = 2.75;
    right_margin = 1.25;
    top_margin = 1;
    bottom_margin = 2.75;
    
    rel_width = 0.15*N_type;
    rel_height = 0.15*N_type;
    
    fig_width = left_margin+right_margin+rel_width;
    fig_height = bottom_margin+top_margin+rel_height;
   
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
        
    sp1 = subplot(1,1,1);
    hold on;
    mu_plot = NaN(N_type, N_type);
    mu_plot(sub2ind([N_type, N_type], r, c)) = mu_diff;
    imagesc(1:N_type, 1:N_type, mu_plot, 'AlphaData',~isnan(mu_plot));
    set(gca, 'YDir', 'normal');
    scatter(c(p<=p_crit), r(p<=p_crit), 10, 'r', 'filled');
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    xtickangle(90);
    set(gca, 'Xtick', 1:N_type, 'xticklabel', cell_type);
    set(gca, 'Ytick', 1:N_type, 'yticklabel', cell_type);    
    title('| \Delta Average S |', 'FontSize', 6);
    cb(1) = colorbar('EastOutside', 'FontSize', 6);
    
    axis tight;
    box off;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin bottom_margin rel_width rel_height]);
    
    cb(1).Position(3) = 0.03;
    
end
    
function [mu, se, df] = get_t_prop(S, N_neighbors, method)

    mu = nanmean(S);
    sd = sqrt(nansum((S-mu).^2)/sum(~isnan(S)));
    
    if (isequal(method, 'None'))
        df = sum(~isnan(S))-1;
    elseif (isequal(method, 'Global'))
        df = nanmedian(sum(~isnan(S))./N_neighbors)-1;
    elseif (isequal(method, 'Local'))
        df = nanmedian(sum(~isnan(S))./N_neighbors)-1;
    end
    
    se = sd/sqrt(df+1);

end

function [p, mu_diff, t_stat, df_joint] = welch_t_test(mu1, mu2, se1, se2, df1, df2)

    mu_diff = abs(mu1-mu2);
    se_del = sqrt(se1^2+se2^2);
    t_stat = mu_diff/se_del;
    df_joint = se_del^4/(se1^4/df1+se2^4/df2);
    p = 2*(1-tcdf(t_stat, df_joint));
    
    if (floor(df1) <= 5 || floor(df2) <= 5)
        t_stat = NaN;
        p = NaN;
    end
    
end