function write_mtx(CBs, genes, counts, filepath)

    if (exist(filepath, 'dir')~=7)
        mkdir(filepath);
    end 

    cb_file = sprintf('%s/barcodes.tsv', filepath);
    fid = fopen(cb_file, 'wt');
    fprintf(fid, '%s\n', CBs{:});
    fclose(fid);
    gzip(cb_file);
    delete(cb_file);
    
    gene_file = sprintf('%s/features.tsv', filepath);
    fid = fopen(gene_file, 'wt');
    for i = 1:length(genes)
        fprintf(fid, '%s\t%s\t%s\n', genes{i}, genes{i}, 'Gene Expression');
    end
    fclose(fid);
    gzip(gene_file);
    delete(gene_file);
    
    [r, c, v] = find(counts);
    mtx_file = sprintf('%s/matrix.mtx', filepath);
    fid = fopen(mtx_file, 'wt');    
    
    fprintf(fid, '%s\n', '%%MatrixMarket matrix coordinate integer general');
    fprintf(fid, '%s\n', '%metadata_json: {"format_version": 2, "software_version": "3.1.0"}');
    
    fprintf(fid, '%d %d %d\n', length(genes), length(CBs), length(v));
    for i = 1:length(v)
        fprintf(fid, '%d %d %d\n', c(i), r(i), v(i));
    end
    fclose(fid);
    gzip(mtx_file);
    delete(mtx_file);
    
end