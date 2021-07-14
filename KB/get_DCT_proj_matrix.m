function [C, alpha] = get_DCT_proj_matrix(include_const, max_t_fourier, max_p_fourier, max_t_pow, max_p_pow)

    % A. Constant terms
    
    if (include_const==true)
        C.const.N = 8;
        C.const.mat = eye(C.const.N);
    else
        C.const.N = 0;
        C.const.mat = zeros(8,0);
    end

    % B. Terms that are periodic in only theta or phi, that satisfy
    % boundary conditions
    
    C.fourier.max_t = max_t_fourier;
    C.fourier.theta = 2:2:max_t_fourier;
    C.fourier.N_theta = length(C.fourier.theta);
    
    C.fourier.max_p = max_p_fourier;
    C.fourier.phi   = 1:max_p_fourier;
    C.fourier.N_phi   = length(C.fourier.phi);    
    
    C.fourier.mat = repmat(eye(8), 1, 2*C.fourier.N_theta+C.fourier.N_phi);

    % C. Terms that are periodic in both theta or phi, satisfying boundary
    % conditions
   
    valid_ind = fliplr(triu(ones(max_t_pow+1)));
    valid_ind(1,1) = 0;
    
    [order_y, order_x] = find(valid_ind);
    order_y = order_y-1;
    order_x = order_x-1;
    valid_ind = find(valid_ind);
    poly_parity = mod(order_y+order_x,2);
    
    vals = sortrows([poly_parity, order_y+order_x, order_y, order_x, valid_ind], ...
                    {'ascend'; 'ascend'; 'descend'; 'ascend'; 'ascend'});
   
    odd_poly.tot_order = vals(vals(:,1)==1,2);
    assert(all(mod(odd_poly.tot_order,2)==1));
    odd_poly.order_y   = vals(vals(:,1)==1,3);
    odd_poly.order_x   = vals(vals(:,1)==1,4);
    odd_poly.mat_ind   = vals(vals(:,1)==1,5);    
     
    even_poly.tot_order = vals(vals(:,1)==0,2);
    assert(all(mod(even_poly.tot_order,2)==0));
    even_poly.order_y   = vals(vals(:,1)==0,3);
    even_poly.order_x   = vals(vals(:,1)==0,4);
    even_poly.mat_ind   = vals(vals(:,1)==0,5);
            
    odd_poly.contrib = [1,   0, 2/3,   0;
                        0,   1,   0, 2/3;
                        0,   0,   0,   0;
                        0,   0,   0,   0;
                        0,   0,   0,   0;
                        0,   0,   1,   0;
                        0,   0,   0,   1;
                        0,   0,   0,   0];
                    
    odd_poly.parity_type = zeros(length(odd_poly.mat_ind),1);
                    
    odd_poly.parity_type(mod(odd_poly.order_y,2)==1 & odd_poly.order_x==0) = 1;                                % odd_only_y
    odd_poly.parity_type(mod(odd_poly.order_x,2)==1 & odd_poly.order_y==0) = 2;                                % odd_only_x
    odd_poly.parity_type(mod(odd_poly.order_y,2)==1 & odd_poly.order_x >0 & mod(odd_poly.order_x,2)==0 ) = 3;  % odd_odd_y_even_x
    odd_poly.parity_type(mod(odd_poly.order_y,2)==0 & odd_poly.order_y >0 & mod(odd_poly.order_x,2)==1 ) = 4;  % odd_odd_x_even_y
               
    even_poly.contrib = [0,   0,   0,   0;
                         0,   0,   0,   0;
                         1,   0,   0, 2/3;
                         0,   1,   0, 2/3;
                         0,   0,   1,   0;
                         0,   0,   0,   0;
                         0,   0,   0,   0;
                         0,   0,   0,   1];
                
    even_poly.parity_type = zeros(length(even_poly.mat_ind),1);
    
    even_poly.parity_type(mod(even_poly.order_y,2)==0 & even_poly.order_x==0) = 1;                                % even_only_y
    even_poly.parity_type(mod(even_poly.order_x,2)==0 & even_poly.order_y==0) = 2;                                % even_only_x    
    even_poly.parity_type(mod(even_poly.order_y,2)==1 & mod(even_poly.order_x,2)==1) = 3;                         % even_odd_y_odd_x        
    even_poly.parity_type(mod(even_poly.order_y,2)==0 & mod(even_poly.order_x,2)==0 & even_poly.order_y>0 & even_poly.order_x>0) = 4; % even_even_y_even_x
           
    C.poly.coeffs = [ sqrt(6);
                     -sqrt(6);
                      sqrt(6);
                      sqrt(6);
                     -sqrt(8);
                      4*sqrt(3)/3;
                     -4*sqrt(3)/3;
                      2*sqrt(6)/3];
    
    C.poly.odd = odd_poly;
    C.poly.even = even_poly;
    C.poly.t_pow = max_t_pow;
    C.poly.p_pow = max_p_pow;
    C.poly.p_pow_sin = [1:2:C.poly.p_pow];
    C.poly.p_pow_cos = [1:C.poly.p_pow];    
    
    C.poly.mat = C.poly.coeffs.*[repmat(C.poly.odd.contrib (:,C.poly.odd.parity_type),  1, length(C.poly.p_pow_sin)), ...
                                 repmat(C.poly.even.contrib(:,C.poly.even.parity_type), 1, length(C.poly.p_pow_cos))];                  
        
    C.N.const     = C.const.N;
    C.N.t_fourier = 2*C.fourier.N_theta*8;
    C.N.p_fourier =   C.fourier.N_phi*8;
    C.N.poly      = length(C.poly.odd.mat_ind)*length(C.poly.p_pow_sin)+length(C.poly.even.mat_ind)*length(C.poly.p_pow_cos);
    C.N.total     = C.N.const+C.N.t_fourier+C.N.p_fourier+C.N.poly;
    
    alpha = zeros(1, C.N.total);
    
    offset = C.N.const+C.N.t_fourier+C.N.p_fourier;    
    C.init_alpha_ind.y  = offset+find(C.poly.odd.order_y==1 & C.poly.odd.order_x==0);
    C.init_alpha_ind.x  = offset+find(C.poly.odd.order_x==1 & C.poly.odd.order_y==0);
    
    offset = offset+length(C.poly.odd.mat_ind)*length(C.poly.p_pow_sin);
    
    C.init_alpha_ind.xy = offset+find(C.poly.even.order_y==1 & C.poly.even.order_x==1);
    C.init_alpha_ind.y2 = offset+find(C.poly.even.order_y==2 & C.poly.even.order_x==0);
    C.init_alpha_ind.x2 = offset+find(C.poly.even.order_y==0 & C.poly.even.order_x==2);
    
    C.init_alpha_val.y  = 1;
    C.init_alpha_val.x  = 1;
    C.init_alpha_val.xy = 2;
    C.init_alpha_val.y2 = 1;
    C.init_alpha_val.x2 = 1;
    
    alpha(struct2array(C.init_alpha_ind)) = struct2array(C.init_alpha_val);
    
end

