function paper_print(filename, opts)

    if (nargin == 1)        
        print('-dtiff','-r600', sprintf('%s.tiff', filename));        
    elseif (strcmp(opts, 'SVG'))
        print(sprintf('%s.svg', filename), '-dsvg');
    elseif (strcmp(opts, 'PNG'))
        print(sprintf('%s.png', filename), '-dpng', '-r2400');
    end
    
end
    
        
