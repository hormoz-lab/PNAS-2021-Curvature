function [X, vars] = rejection_sample(X, N_points, targ_density, bounds, skew)

    assert(isequal(size(bounds), [2, 2]));   
    if (nargin == 4)
        skew = 10;
    end

    N_rejection = N_points*skew;
    
    var1_sweep = bounds(1,1)+(bounds(1,2)-bounds(1,1))*rand(N_rejection,1);
    var2_sweep = bounds(2,1)+(bounds(2,2)-bounds(2,1))*rand(N_rejection,1);
    
    targ_sweep = arrayfun(@(v1,v2) targ_density(v1,v2), var1_sweep, var2_sweep);
    targ_trial = max(targ_sweep)*rand(N_rejection,1);
    
    mask = targ_trial <= targ_sweep;
    var1_sweep = var1_sweep(mask);
    var2_sweep = var2_sweep(mask);
    
    assert(length(var1_sweep) > N_points);
    
    var1_sweep = var1_sweep(1:N_points);
    var2_sweep = var2_sweep(1:N_points);
    
    X = matlabFunction(X, 'Vars', symvar(X));
    
    X = arrayfun(@(v1,v2) X(v1,v2), var1_sweep, var2_sweep, 'un', false);
    X = [horzcat(X{:})]';
    
    vars = [var1_sweep var2_sweep];    
end


