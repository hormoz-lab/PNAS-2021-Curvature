function preprocess_Hateren(data_path, results_dir, do_full)

    % do_full = true needs a bigmem machine, the final saved file is ~130GB

    N_images = 4212;

    hc_patches = cell(N_images, 1);
    hc_patch_index = cell(N_images, 1);
    DCT = load(sprintf('%s/KB/DCTParams.mat', results_dir));

    width = 1536; height = 1024; offset = 2;
    [centre_c, centre_r] = meshgrid(2:width-2*offset-1, 2:height-2*offset-1);
    if (do_full)
        N_subsample = numel(centre_c);
    else
        N_subsample = 5000;
    end
    contrast_pct = 0.20;
    rng(3439);

    for i = 1:N_images
        fprintf('Processing image %d of %d\n', i, N_images);        
        f1 = fopen(sprintf('%s/Images/imk%05d.iml', data_path, i), 'rb', 'ieee-be');
        if (f1 == -1)
            continue;
        end
        buf = fread(f1, [width, height], 'uint16')';
        fclose(f1);    
        buf = log(buf(offset+1:height-offset, offset+1:width-offset)+1);    
        assert(all(isfinite(buf(:))), 'Zero pixel in image %d', i);
        
        if (do_full)
            subsample_ind = [1:N_subsample]';
        else
            subsample_ind = randperm((width-3*offset)*(height-3*offset), N_subsample);
        end
        patches = arrayfun(@(r,c) buf(r-1:r+1, c-1:c+1), centre_r(subsample_ind), centre_c(subsample_ind), 'un', false);
        patches = cellfun(@(x) reshape(x, 9, 1), patches, 'un', false);
        patches = horzcat(patches{:});

        patches_Dnorm = sqrt(sum((DCT.D*patches).*patches,1));
        [~, ind] = sort(patches_Dnorm, 'descend');
        ind = ind(1:(N_subsample*contrast_pct));
        hc_patch_index{i} = subsample_ind(ind);
        hc_patches{i} = patches(:, ind);    

    end
    
    if (~exist(results_dir, 'dir'))
        mkdir(results_dir);
    end
    
    if (do_full)
        outfile = 'HaterenAll';
    else
        outfile = 'Hateren';
    end
    
    save(sprintf('%s/KB/%s.mat', results_dir, outfile), 'N_images', 'hc_patches', 'hc_patch_index', 'width', 'height', 'offset', ...
        'centre_c', 'centre_r', 'N_subsample', 'contrast_pct', '-v7.3');
    
    clearvars -except hc_patches DCT results_dir outfile
    
    R9 = horzcat(hc_patches{:});
    clear hc_patches;

    Dnorm = sqrt(sum((DCT.D*R9).*R9,1));

    S7 = R9./Dnorm;
    R9 = R9';
    save(sprintf('%s/KB/%s.mat', results_dir, outfile), 'R9', '-append');
    clear R9;

    S7 = DCT.lambda*DCT.A'*(S7);
    S7 = S7';

    save(sprintf('%s/KB/%s.mat', results_dir, outfile), 'S7', '-append');
    clearvars -except S7 results_dir outfile

   if (do_full)
       % This needs to be parallelized to run in a reasonable amount of
       % time and probably needs a bigmem machine - consider this bit pseudo-code
       [~, D_knn_600] = knnsearch(S7, S7, 'K', 600);
       D_knn_600 = D_knn_600(:,end);
       S7 = S7(D_knn_600<=prctile(D_knn_600,10),:);
       [S7_unique, uniq_ind] = unique(S7, 'rows', 'stable');
       save(sprintf('%s/KB/%s.mat', results_dir, outfile), 'D_knn_600', 'S7_unique', 'uniq_ind', '-append');
   else
       [~, D_knn_100] = knnsearch(S7, S7, 'K', 100);
       D_knn_100 = D_knn_100(:,end);
       save(sprintf('%s/KB/%s.mat', results_dir, outfile), 'D_knn_100', '-append');
   end
   
end