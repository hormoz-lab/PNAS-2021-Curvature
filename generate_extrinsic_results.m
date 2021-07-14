function generate_extrinsic_results(results_dir)
        
    N = 10000;
 
    %% Sample S^n    
    
    rng(72839); % Same seed as for intrinsic    
    n = 2;
    data = ManifoldHandler.sample_hypersphere(N,n);   
    manifold_curvature(sprintf('%s/ToyData/S%d', results_dir, n), data, n, 0.017);
    
    rng(9829);
    n = 3;
    data = ManifoldHandler.sample_hypersphere(N,n);    
    manifold_curvature(sprintf('%s/ToyData/S%d', results_dir, n), data, n, 0.02);
    
    rng(34124);
    n = 5;
    data = ManifoldHandler.sample_hypersphere(N,n);    
    manifold_curvature(sprintf('%s/ToyData/S%d', results_dir, n), data, n, 0.028);
    
    rng(58305049);
    n = 7;
    data = ManifoldHandler.sample_hypersphere(N,n);    
    manifold_curvature(sprintf('%s/ToyData/S%d', results_dir, n), data, n, 0.055);                
     
    %% Sample H2
    
    rng(5849989);
    % Pick (a,b,c) so S_min = -2
    a = 2; b = 2; c = 1;
    data = ManifoldHandler.sample_one_sheet_hyperboloid(3.25*N, a, b, c);
    mask = abs(data(:,3))<1;
    assert(sum(mask)>=N);
    in_mask = find(mask, N, 'first');
    out_mask = find(~mask, round(N/sum(mask)*(sum(~mask))), 'first');
    data = [sortrows(data(in_mask,:), 3); data(out_mask,:)];
    analytic = ManifoldHandler.one_sheet_hyperboloid_curvature(data(1:N,3), a, b, c);
    manifold_curvature(sprintf('%s/ToyData/H2', results_dir), data, 2, 0.022, 'calib_points', [1:N]', 'point_ids', [1:N]');
    save(sprintf('%s/ToyData/H2/Curvature.mat', results_dir), 'a', 'b', 'c', 'analytic', '-append');
        
     %% Sample T2
     
    rng(5849);
    % Pick (R,r) so that S_min = -2 => S_max = 4/3;
    R = 5/2;
    r = 1/2;
    [data, theta] = ManifoldHandler.sample_torus(N, R, r);
    [~, ind] = sort(theta(:,1));
    data = data(ind,:);
    theta = theta(ind,:);
    analytic = ManifoldHandler.torus_curvature(R, r, theta(:,1));
    manifold_curvature(sprintf('%s/ToyData/T2', results_dir), data, 2, 0.03);
    save(sprintf('%s/ToyData/T2/Curvature.mat', results_dir), 'theta', 'R', 'r', 'analytic', '-append');
    
end