function W = gaussian_weight_matrix(data, kernel_eps)

    N = size(data,1);
    W = exp(-squareform(pdist(data, 'squaredeuclidean'))/kernel_eps)-eye(N);
    
end

