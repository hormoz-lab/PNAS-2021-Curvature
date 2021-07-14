function idx = knnsearch_block(data, K)
    
    N = size(data,1);        
    N_blocks = ceil(N*K/1e8);
    block_size = ceil(N/N_blocks);
    
    idx = cell(N_blocks,1);
    
    for i = 1:N_blocks
        [i N_blocks]
        lo = (i-1)*block_size+1;
        hi = min(i*block_size, N);
        if (i == N_blocks)
            assert(hi == N);
        end        
        idx{i} = knnsearch(data, data(lo:hi,:), 'K', K);        
    end
    
    idx = vertcat(idx{:});    
end