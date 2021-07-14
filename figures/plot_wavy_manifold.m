function plot_wavy_manifold()

    fig_width = 5.4;
    fig_height = 3.8;
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    [X, Y, Z] = peaks;
    X = X-min(X(:));
    Y = Y-min(Y(:));
    Z = Z/max(abs(Z(:)));
    L = max(sqrt(X(:).^2+Y(:).^2))+0.1;
    Y = Y./L;
    X = X./L;
    Z = Z+sqrt(1-X.^2-Y.^2);

    rng(0890);
    sig = 0.01;
    Xp = X.*Z + sig*randn(size(X));
    Yp = X.*Y + sig*randn(size(X));
    Zp = Y.*Y + sig*randn(size(X));

    mask = Xp>0.55 & Zp>0.2;
    ind = find(Xp == max(max(Xp)));
    del = 0.2;    
           
    sp1 = subplot(1,1,1);
    surf(Z.*X, X.*Y, Y.*Y, 0.2*ones([size(X), 3]), 'FaceAlpha', 0.1, 'Edgecolor', [0.6 0.6 0.6]);
    hold on;
    scatter3(Xp(:), Yp(:), Zp(:), 2, 'b', 'filled'); colorbar;
    scatter3(Xp(mask), Yp(mask), Zp(mask), 3, 'g', 'filled');    
    fill3(Xp(ind)*ones(1,4), Yp(ind)+[0 del del 0], Zp(ind)+[0 0 del del], 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'none');    
    nudge = 0.04;
    line(Xp(ind)+[0 del], Yp(ind)*ones(1,2), Zp(ind)*ones(1,2), 'Color', 'c', 'LineWidth', 1.0);
    scatter3(Xp(ind), Yp(ind), Zp(ind), 9, 'r', 'filled');
    text(Xp(ind), Yp(ind), Zp(ind)-nudge, '$p$', 'FontSize', 5, 'interpreter', 'latex');
    text(Xp(ind)-nudge, Yp(ind)+del/2-nudge, Zp(ind)+del/2, '$T_M(p)$', 'FontSize', 5, 'interpreter', 'latex');
    text(Xp(ind)+del+nudge, Yp(ind), Zp(ind), '$N_M(p)$', 'FontSize', 5, 'interpreter', 'latex');
    hold off;    
    view([30 30]);
    
    grid off;
    axis tight;
    box off;    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'ZAxis'), 'FontSize', 5);
    set(gca, 'LineWidth', 1.0, 'TickDir', 'Out');
    xlabel('x', 'FontSize', 6);
    ylabel('y', 'FontSize', 6);
    zlabel('z', 'FontSize', 6);    

    left_margin = 0.7;
    right_margin = 0.3;    
    top_margin = 0.1;
    bottom_margin = 0.5;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);

end
