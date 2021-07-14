function plot_KB_patches(Klein_fcn, x_label, y_label)

    fig_width = 6.6;
    fig_height = 6.6;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
      
    theta_grid = linspace(0,pi,9);
    N_theta = length(theta_grid);
    phi_grid = linspace(0,2*pi,9);
    N_phi = length(phi_grid);
    [theta_sweep, phi_sweep] = meshgrid(theta_grid, phi_grid);
    theta_sweep = theta_sweep(:);
    phi_sweep = phi_sweep(:);
    [~, where_phi] = ismember(phi_sweep, phi_grid);
    [~, where_theta] = ismember(theta_sweep, theta_grid);
        
    patches = arrayfun(@(t,p) Klein_fcn(t, p)', theta_sweep, phi_sweep, 'un', false);
    patches = cellfun(@(x) (x-min(x)), patches, 'un', false);
    patches = cellfun(@(x) reshape(x/max(x), [3, 3]), patches, 'un', false);
    enlarge_fac = 5;
    half_width = ceil(3*enlarge_fac/2);
    offset = 10;
    patch_grid = ones(2+N_phi*3*enlarge_fac+(N_phi-1)*offset, 2+N_theta*3*enlarge_fac+(N_theta-1)*offset,3);
    
    patches = cellfun(@(x) repelem(x, enlarge_fac, enlarge_fac), patches, 'un', false);
    
    border_color = zeros(1,1,3);
    border_color(1,1,:) = [0.5 0 0.5];
    
    for i = 1:length(phi_sweep)
        theta_centre = 1+half_width+(3*enlarge_fac+offset)*(where_theta(i)-1);
        phi_centre   = 1+half_width+(3*enlarge_fac+offset)*(where_phi(i)-1);        
        patch_grid(phi_centre+[-half_width  :half_width  ],theta_centre+[-half_width  :half_width  ],:) = repmat(border_color, [2*half_width+1 2*half_width+1, 1]);
        patch_grid(phi_centre+[-half_width+1:half_width-1],theta_centre+[-half_width+1:half_width-1],:) = repmat(patches{i}, [1,1,3]);
    end

    sp1 = subplot(1,1,1);
    hold on;    
    h = imagesc(patch_grid);
    colormap gray;
    
    a_len = 2;
    a_wid = 2;
    
    for i = 1:(N_phi-1)
        y_start = 1.5+enlarge_fac*3+(i-1)*(offset+enlarge_fac*3);
        y_end   = 0.5+(enlarge_fac*3+offset)*i;        
        line([1+half_width 1+half_width], [y_start+1 y_end], 'Color', 'r', 'LineWidth', 1.0);
        line([size(patch_grid,2)-half_width size(patch_grid,2)-half_width], [y_start+1 y_end], 'Color', 'r', 'LineWidth', 1.0);
        if (mod(i-1,3)==0)
            patch([1+half_width-a_wid 1+half_width 1+half_width+a_wid], 1+[y_start+a_len y_end-a_len y_start+a_len], 'r', 'EdgeColor', 'r');
            patch([size(patch_grid,2)-half_width-a_wid size(patch_grid,2)-half_width size(patch_grid,2)-half_width+a_wid], [y_end-a_len y_start+a_len y_end-a_len], 'r', 'EdgeColor', 'r');
        end
    end
    
    a_wid = 2;
    for i = 1:(N_theta-1)
        x_start = 1.5+enlarge_fac*3+(i-1)*(offset+enlarge_fac*3);
        x_end   = 0.5+(enlarge_fac*3+offset)*i;        
        line([x_start+1 x_end], [1+half_width 1+half_width], 'Color', 'b', 'LineWidth', 1.0);
        line([x_start+1 x_end], [size(patch_grid,1)-half_width size(patch_grid,1)-half_width], 'Color', 'b', 'LineWidth', 1.0);
        if (mod(i-1,3)==0)
            patch(1+[x_start+a_len x_end-a_len x_start+a_len], [1+half_width-a_wid 1+half_width 1+half_width+a_wid], 'b', 'EdgeColor', 'b');
            patch(1+[x_start+a_len x_end-a_len x_start+a_len], [size(patch_grid,2)-half_width-a_wid size(patch_grid,2)-half_width size(patch_grid,2)-half_width+a_wid], 'b', 'EdgeColor', 'b');
        end
    end
        
    hold off;
    
    axis tight;    
    xticks([]);
    yticks([]);
    box off;
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'ZAxis'), 'FontSize', 5);
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    
    xlabel(x_label, 'FontSize', 6, 'interpreter', 'latex');
    ylabel(y_label, 'FontSize', 6, 'interpreter', 'latex');    
    
    ax = gca;
    axis(ax,'off')    
    ax.XLabel.Visible = 'on';
    ax.YLabel.Visible = 'on';
    
    left_margin = 0.4;
    right_margin = 0.1;    
    top_margin = 0.1;
    bottom_margin = 0.4;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);

end
