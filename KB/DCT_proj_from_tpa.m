function out = DCT_proj_from_tpa(C, theta, phi, alpha)

    % A. Constant terms   
    const_terms = [];
    if (C.const.N > 0)
        const_terms = ones(1,C.const.N);        
    end
    
    % B. Terms that are periodic in only theta or phi, that satisfy
    % boundary conditions
    
    fourier_terms = [];
    if (C.N.t_fourier > 0)
        fourier_terms = repelem([cos(theta*C.fourier.theta) sin(theta*C.fourier.theta)], 1, 8);
    end
    if (C.N.p_fourier > 0)
        fourier_terms = [fourier_terms repelem(cos(phi*C.fourier.phi), 1, 8)];
    end

    % C. Terms that are periodic in both theta or phi, satisfying boundary
    % conditions
    
    a = cos(theta);
    b = sin(theta);
    c = cos(phi);
    d = sin(phi);
    
    theta_matrix = (b.^[0:C.poly.t_pow])'*(a.^[0:C.poly.t_pow]);
    
    poly_odd_terms  = theta_matrix(C.poly.odd.mat_ind) .*(d.^[C.poly.p_pow_sin]);
    poly_even_terms = theta_matrix(C.poly.even.mat_ind).*(c.^[C.poly.p_pow_cos]);
        
    poly_terms = [poly_odd_terms(:); poly_even_terms(:)]';    
              
    out = [C.const.mat C.fourier.mat C.poly.mat].*[const_terms fourier_terms poly_terms];        
    assert(size(out,2) == C.N.total);
    if (nargin == 4)
        assert(length(alpha)==C.N.total);
        out = alpha*out';
    end
end