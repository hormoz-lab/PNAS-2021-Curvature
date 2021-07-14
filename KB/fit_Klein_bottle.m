function [alpha, v_est, v_err, X] = fit_Klein_bottle(v, C, theta, phi, v_norm)

    L = size(v,1);
    assert(length(theta) == L);
    assert(length(phi)   == L);
                   
    y = v';
    y = y(:);
    
    X = arrayfun(@(t,p) DCT_proj_from_tpa(C, t, p), theta, phi, 'un', false);
    X = vertcat(X{:});
  
    abs_tol = 1e-2;
    rel_tol = 1e-3;
    N_trials = 100;
    i = 1;
    v_err = 1e10;
    v_err_prev = inf;
    v_est = [];
    alpha = 0;
            
    if (nargin == 4)        
        assert(length(v_norm)== L);
        X = X./repelem(v_norm, 8);
    end    
    
    while ((i <= N_trials) && (sum(v_err) > abs_tol) && (sum(v_err_prev)-sum(v_err))>rel_tol*sum(v_err))
        
        XtX = X'*X;
        Xty = X'*y;

        alpha_prev = alpha;
        alpha = (XtX)\Xty;

        v_est = X*alpha;
        v_est = reshape(v_est, [8, L])';
        v_norm = sqrt(sum(v_est.^2,2));
        v_est = v_est./v_norm;
        v_err_prev = v_err;
        v_err = sqrt(sum((v-v_est).^2,2));
        
        if (sum(v_err_prev) < sum(v_err))
            alpha = alpha_prev;
            v_err = v_err_prev;
            break;
        else        
            fprintf('%d: Err=%g, %s\n', i, sum(v_err), num2str(alpha'));
            X = X./repelem(v_norm, 8);
            i = i+1;
        end
    end      
end
