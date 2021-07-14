function subsample_transcripts(data_path)

    [counts, CBs, genes] = load_10x_mtx(data_path);
    
    [r, c, v] = find(counts);
    
    rng(78242483);
    data = repelem([1:length(r)]', v);
    length(data)
    
    N_ss = 2;
    
    for i = 1:N_ss
        fprintf('Downsampling %d\n', i);
        if (i==1)
            ss{i} = datasample(data, round(length(data)/2), 'Replace', false);
        else
            ss{i} = datasample(ss{i-1}, round(length(ss{i-1})/2), 'Replace', false);
        end
        length(ss{i})        
    end    
    
    for i = 1:N_ss
        fprintf('Matrix %d\n', i);
        count{i} = accumarray([r c], accumarray(ss{i}, 1, size(r)), [length(CBs) length(genes)], [], [], true);      
        if (i==1)
            assert(isequal(logical(counts(:)) & logical(count{i}(:)), logical(count{i}(:))));            
        else
            assert(isequal(logical(count{i-1}(:)) & logical(count{i}(:)), logical(count{i}(:))));                        
        end
        mkdir(sprintf('%s/SSTranscripts_Fact%d', data_path, i));
        write_mtx(CBs, genes, count{i}, sprintf('%s/SSTranscripts_Fact%d', data_path, i));
    end
    
end