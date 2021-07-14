function [theta_est, phi_est, norm_est] = get_Klein_coords_opt(v, C, alpha, theta_init, phi_init)
    
    tic; 
    
    assert(size(v,2)==8);
    assert(length(alpha)==C.N.total);
    L = size(v,1);
    
    if (nargin > 3)
        assert(nargin == 5);
    else    
        [theta_init, phi_init] = get_Klein_coords_init(v, C, alpha);
    end
    
    f = @(x,v) err_func(x,v,C,alpha);
    
    sol = zeros(L,2);
    resnorm = zeros(L,1);
    warning('off', 'optim:fsolve:NonSquareSystem');
    warning('off', 'optimlib:checkbounds:IgnoringExtraLbs');
    options = optimoptions(@lsqnonlin,'Display','off', 'StepTolerance', 1e-3, 'FunctionTolerance', 1e-6, 'UseParallel', false);
    
    parfor i = 1:L
        fprintf('Optimized %d/%d\n', i, L);
        [sol(i,:), resnorm(i)] = lsqnonlin(@(x) f(x,v(i,:)), [theta_init(i); phi_init(i)], ...
                                                [0; 0], [pi; 2*pi], options);
    end
        
    theta_est = sol(:,1);
    phi_est   = sol(:,2);
    
    fprintf('(Theta/phi) optimization ran in %g seconds. Err=%g\n', toc, sum(sqrt(resnorm)));
    
    out = arrayfun(@(t,p) DCT_proj_from_tpa(C, t, p, alpha), theta_est, phi_est, 'un', false);
    out = vertcat(out{:});
    norm_est = sqrt(sum(out.^2,2));   
    
end

function f = err_func(x, v, C, alpha)
    f = DCT_proj_from_tpa(C, x(1), x(2), alpha);
    f = f/norm(f)-v;
end