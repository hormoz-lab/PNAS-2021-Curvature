function make_intro_subplots(results_dir)

    outdir = sprintf('%s/Figures/Intro', results_dir);
    if (~exist(outdir, 'dir'))
        mkdir(outdir);
    end

    load(sprintf('%s/Intrinsic.mat', results_dir), 'X');
    
    plot_manifold_with_S(X, X(:,3), [min(X(:,3)), max(X(:,3))], '10K Points from $\mathcal{S}^2$', 'z');
    paper_print(sprintf('%s/Sphere', outdir));
    
    plot_flattened_manifold([atan2(X(:,2), X(:,1))+pi asin(X(:,3))], X(:,3), 'Polar Coordinates', '\theta_1 (rad)', '\theta_2 (rad)', ...
                            [0, 2*pi], [-pi/2, pi/2], [0, 2*pi], [-pi/2, pi/2], {'0'; '2$\pi$'}, {'$-\frac{\pi}{2}$'; '$\frac{\pi}{2}$'});
    paper_print(sprintf('%s/Sphere_Polar', outdir));
    
    P = 1.0./(1-X(:,3)).*[X(:,1:2)];
    mask = abs(P(:,1))<5 & abs(P(:,2))<5;
    plot_flattened_manifold(P(mask,:), X(mask,3), 'Stereographic Projection', 'x', 'y', [-5 5], [-5 5], [-1 1], [-1 1]);
    paper_print(sprintf('%s/Sphere_Stereographic', outdir));
    
    N = 10000;
    rng(394893);
    P = [rand(N,1)*2 rand(N,1)*20];
    
    a = 1;
    b = 1/(2*pi);
    
    opts = optimoptions('fsolve','Display','off');
    theta = arrayfun(@(L) fsolve(@(tmax) integral(@(t) sqrt(a^2+b^2+2*a*b*t+b^2*t.^2), 0, tmax)-L, 0, opts), P(:,2));
    R = a+b*theta;
    Y = R.*cos(theta);
    Z = R.*sin(theta);
    
    plot_manifold_with_S([P(:,1) Y Z], P(:,2), [min(P(:,2)), max(P(:,2))], '10K Points from Scroll', 'Arclength');
    paper_print(sprintf('%s/Scroll', outdir));
        
    plot_flattened_manifold(P, P(:,2), 'Unfurled Scroll', 'x', 'Arclength', [0 2], [0 20], [0 2], [0 20]);
    paper_print(sprintf('%s/Scroll_Unfurled', outdir));

    close all;
    
end