function [X, vars] = uniform_sample(X, N_points, bounds)

    g = ManifoldHandler.get_metric(X);
    det_g_fcn = matlabFunction(simplify(sqrt(det(g))), 'Vars', symvar(X));
    
    [X, vars] = ManifoldHandler.rejection_sample(X, N_points, det_g_fcn, bounds);
    
end


