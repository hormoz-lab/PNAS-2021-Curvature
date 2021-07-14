function setup_DCT_basis(results_dir)

    A = [[ 1  0 -1  1  0 -1  1  0 -1] / sqrt(6)  ;
         [ 1  1  1  0  0  0 -1 -1 -1] / sqrt(6)  ;
         [ 1 -2  1  1 -2  1  1 -2  1] / sqrt(54) ;
         [ 1  1  1 -2 -2 -2  1  1  1] / sqrt(54) ;
         [ 1  0 -1  0  0  0 -1  0  1] / sqrt(8)  ;
         [ 1  0 -1 -2  0  2  1  0 -1] / sqrt(48) ;
         [ 1 -2  1  0  0  0 -1  2 -1] / sqrt(48) ;
         [ 1 -2  1 -2  4 -2  1 -2  1] / sqrt(216)]';

    lambda = diag(1./sum(A.^2,1));   
    
    % This is just D = A*lambda^2*A', but don't want to keep things integer
    
    D = [ 2 -1  0 -1  0  0  0  0  0;
         -1  3 -1  0 -1  0  0  0  0;
          0 -1  2  0  0 -1  0  0  0;
         -1  0  0  3 -1  0 -1  0  0;
          0 -1  0 -1  4 -1  0 -1  0;
          0  0 -1  0 -1  3  0  0 -1;
          0  0  0 -1  0  0  2 -1  0;
          0  0  0  0 -1  0 -1  3 -1;
          0  0  0  0  0 -1  0 -1  2];
    
    if (~exist(results_dir, 'dir'))
        mkdir(results_dir);
    end
    
    save(sprintf('%s/DCTParams.mat', results_dir), 'A', 'D', 'lambda');
    
end
