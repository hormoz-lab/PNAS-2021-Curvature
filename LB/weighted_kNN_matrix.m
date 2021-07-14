function W = weighted_kNN_matrix(data, kernel_eps, k)
    
    N = size(data,1);
    assert(k < N);

    W = gaussian_weight_matrix(data, kernel_eps);
    
    idx = knnsearch(data, data, 'K', k+1);
    [r, ~, c] = find(idx(:,2:end));
    
    mask = logical(full(sparse(r, c, 1, N, N)));
    mask = mask | mask';
    W(~mask) = 0;
    
end

