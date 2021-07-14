function generate_intrinsic_results(results_dir)

    if (~exist(results_dir, 'dir'))
        mkdir(results_dir);
    end
    
    fprintf('Sampling S2...\n');
    
    rng(72839);
    N = 10000;    
    X = ManifoldHandler.sample_hypersphere(N, 2);
    kernel_eps = get_kernel_width(N);

    fprintf('Computed weight matrices...\n');
    
    W_gauss  = gaussian_weight_matrix(X, kernel_eps);
    W_1000   = weighted_kNN_matrix   (X, kernel_eps, 1000);
    W_100    = weighted_kNN_matrix   (X, kernel_eps, 100);    
    W_rgraph = neighborhood_matrix   (X, kernel_eps);

    fprintf('Random walk normalization...\n');

    L_gauss  = random_walk_normalization(W_gauss,  kernel_eps);
    L_1000   = random_walk_normalization(W_1000,   kernel_eps);
    L_100    = random_walk_normalization(W_100,    kernel_eps);    
    L_rgraph = random_walk_normalization(W_rgraph, kernel_eps);
    clear W*

    fprintf('Computing eigenspectra...\n');
    
    E_gauss  = sort(eig(L_gauss));
    E_1000   = sort(eig(L_1000));
    E_100    = sort(eig(L_100));    
    E_rgraph = sort(eig(L_rgraph));

    k = [1:1000]';
    true_E = repelem((k-1).*k, 2*k-1);
    clearvars -except results_dir true_E E_* k
    
    save(sprintf('%s/Intrinsic.mat', results_dir), 'k', 'true_E', 'E*');
    
    fprintf('Heat traces...\n');
    
    n1 = 1000;
    n2 = 49;
    x_lo = 0;
    x_hi = 1.5;
    N_x  = x_hi/0.01;
    x    = linspace(x_lo,x_hi,N_x)';
    f    = 0.21;
    
    traces{1} = arrayfun(@(x) 4*pi*x^2*sum(exp(- true_E*x^2)), x); traces{1}(traces{1}<4*pi) = 4*pi;
    traces{2} = arrayfun(@(x) 4*pi*x^2*sum(exp(- true_E(1:n1)*x^2)), x);
    traces{3} = arrayfun(@(x) 4*pi*x^2*sum(exp(- true_E(1:n2)*x^2)), x);
    traces{4} = arrayfun(@(x) 4*pi*x^2*sum(exp(-E_gauss(1:n1)*x^2)), x);
    traces{5} = arrayfun(@(x) 4*pi*x^2*sum(exp(-E_gauss(1:n2)*x^2)), x);
    traces{6} = arrayfun(@(x) 4*pi*x^2*sum(exp(-(true_E(1:n2)+f*(E_gauss(1:n2)-true_E(1:n2)))*x^2)), x);
    
    frac_error = (true_E(1:n2)-[E_gauss(1:n2) E_1000(1:n2) E_100(1:n2) E_rgraph(1:n2)])./true_E(1:n2);
    frac_error(isinf(frac_error)) = 0;
        
    [b, e] = find(triu(ones(N_x),1));
    x_mid = (x(b)+x(e))/2;
    options = optimoptions(@lsqnonlin,'Display','off', 'StepTolerance', 1e-3, 'FunctionTolerance', 1e-6, 'UseParallel', false);
    
    fprintf('Fitting p2...\n');
    
    sol2 = cellfun(@(v) arrayfun(@(b,e) lsqnonlin(@(S) 4*pi+S*2*pi/3*x(b:e).^2-v(b:e), 2, [], [], options), ...
                                   b, e), traces, 'un', false);
    sol2 = cellfun(@(v) full(sparse(b,e,v,N_x,N_x)), sol2, 'un', false);
    
    fprintf('Fitting p3...\n');
    
    usol3 = cellfun(@(v) arrayfun(@(b,e) lsqnonlin(@(S) 4*pi+S(1)*2*pi/3*x(b:e).^2+S(2)*x(b:e).^3-v(b:e), [2; 0], [], [], options), ...
                                   b, e, 'un', false), traces, 'un', false);
    usol3 = cellfun(@(x) horzcat(x{:})', usol3, 'un', false);
    sgn3 = cellfun(@(S) sign(x_mid.*(S(:,1)*4*pi/3+3*S(:,2).*x_mid)), usol3(1:3), 'un', false);
    root3= cellfun(@(S) -4*pi*S(:,1)/9./S(:,2), usol3(1:3), 'un', false);
    mask3 = cellfun(@(S,r) S>=0 & (r<=x(b) | r>=x(e)), sgn3, root3, 'un', false);
    sol3(1:3) = cellfun(@(v,m) full(sparse(b(m),e(m),v(m,1),N_x,N_x))+full(sparse(b(~m),e(~m),NaN,N_x,N_x)), usol3(1:3), mask3, 'un', false);
    sol3(4:6) = cellfun(@(v) full(sparse(b,e,v(:,1),N_x,N_x)), usol3(4:6), 'un', false);
    
    fprintf('Fitting p4...\n');
    
    usol4 = cellfun(@(v) arrayfun(@(b,e) lsqnonlin(@(S) 4*pi+S(1)*2*pi/3*x(b:e).^2+S(2)*x(b:e).^3+S(3)*x(b:e).^4-v(b:e), [2; 0; 0], [], [], options), ...
                                   b, e, 'un', false), traces, 'un', false);
    usol4 = cellfun(@(x) horzcat(x{:})',usol4, 'un', false);
    sgn4 = cellfun(@(S) sign(x_mid.*(S(:,1)*4*pi/3+3*S(:,2).*x_mid+4*S(:,3).*x_mid.^2)), usol4(1:3), 'un', false);
    det4 = cellfun(@(S) S(:,2).^2*9-(S(:,1).*S(:,3)*64*pi/3), usol4(1:3), 'un', false);
    root4a = cellfun(@(S,d) real((-3*S(:,2)+sqrt(d))./(8*S(:,3))), usol4(1:3), det4, 'un', false);
    root4b = cellfun(@(S,d) real((-3*S(:,2)-sqrt(d))./(8*S(:,3))), usol4(1:3), det4, 'un', false);
    mask4 = cellfun(@(S,d,r1,r2) S>=0 & (d<0 | ((r1<=x(b) | r1>=x(e)) & (r2<=x(b) | r2>=x(e)))), sgn4, det4, root4a, root4b, 'un', false);
    sol4(1:3) = cellfun(@(v,m) full(sparse(b(m),e(m),v(m,1),N_x,N_x))+full(sparse(b(~m),e(~m),NaN,N_x,N_x)), usol4(1:3), mask4, 'un', false);
    sol4(4:6) = cellfun(@(v) full(sparse(b,e,v(:,1),N_x,N_x)), usol4(4:6), 'un', false);    
    
    fprintf('Saving...\n');
    
    save(sprintf('%s/Intrinsic.mat', results_dir));
    
end