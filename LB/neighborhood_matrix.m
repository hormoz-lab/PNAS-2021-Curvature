function W = neighborhood_matrix(data, r)

    N = size(data,1);
    W = (squareform(pdist(data, 'squaredeuclidean'))<r)-eye(N);
    
end

