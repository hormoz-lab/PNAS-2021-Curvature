function plot_K2_in_R3()

    fig_width = 3.5;
    fig_height = 3.4;

    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);

    syms u v;    

    a = @(u) 6*cos(u)*(1+sin(u));
    b = @(u) 16*sin(u);
    c = @(u) 4*(1-cos(u)/2);
    x = piecewise(u <= pi, a(u)+c(u)*cos(u)*cos(v), ...
                  u >  pi, a(u)+c(u)*cos(v+pi));
    y = piecewise(u <= pi, b(u)+c(u)*sin(u)*cos(v), ...
                  u >  pi, b(u));
    z = c(u)*sin(v);          

    x0 = arrayfun(@(t,p) double(subs(x, [u v], [t p])), linspace(0,2*pi, 1000)', 3*pi/2*ones(1000,1));
    y0 = arrayfun(@(t,p) double(subs(y, [u v], [t p])), linspace(0,2*pi, 1000)', 3*pi/2*ones(1000,1));
    z0 = arrayfun(@(t,p) double(subs(z, [u v], [t p])), linspace(0,2*pi, 1000)', 3*pi/2*ones(1000,1));

    colormap('spring');
    sp1 = subplot(1,1,1);
    h = fsurf(x,y,z, [0 2*pi 0 2*pi], 'LineWidth', 0.1);
    hold on;
    line(x0, y0, z0, 'Color', 'b', 'LineWidth', 5);
    hold off;

    camlight(110,70)
    brighten(0.8)
    h.FaceAlpha = 0.75;
    h.AmbientStrength = 0.4;
    view(25, -55);

    xticks([]);
    yticks([]);    
    box off;
        
    title('$K^2 \subset R^3$', 'FontSize', 6', 'FontWeight', 'normal', 'interpreter', 'latex');

    axis tight;
    ax = gca;
    axis(ax,'off')
    
    left_margin = -0.5;
    right_margin = -0.2;    
    top_margin = 0.3;
    bottom_margin = -0.5;
    
    set(sp1, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);

end
