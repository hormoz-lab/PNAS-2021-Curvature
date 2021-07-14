function L = random_walk_normalization(W, kernel_eps)

    N = size(W,1);
    assert(size(W,2)==N);    
    L = 4*(eye(N)-bsxfun(@rdivide, W, sum(W,2)))/kernel_eps;
    
end

