function g = get_metric(X)

    assert(isa(X, 'sym'));
    vars = symvar(X);
    assert(length(vars)==2);

    dX_d1 = diff(X, vars(1));
    dX_d2 = diff(X, vars(2));

    g11 = simplify(sum(dX_d1.*dX_d1));
    g12 = simplify(sum(dX_d1.*dX_d2));
    g21 = simplify(sum(dX_d2.*dX_d1));
    g22 = simplify(sum(dX_d2.*dX_d2));
    
    g = [g11, g12; g21, g22];    
end