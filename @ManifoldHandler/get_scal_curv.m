function S = get_scal_curv(X)

    assert(isa(X, 'sym'));
    vars = symvar(X);
    assert(length(vars)==2);

    g = ManifoldHandler.get_metric(X);

    D1 = simplify(cat(3, diff(g, vars(1)), diff(g, vars(2))));
    ig = simplify(inv(g));    
    C = cat(3, [ManifoldHandler.Christoffel(ig,D1,1,1,1), ManifoldHandler.Christoffel(ig,D1,1,2,1); 
                ManifoldHandler.Christoffel(ig,D1,2,1,1), ManifoldHandler.Christoffel(ig,D1,2,2,1)], ...
               [ManifoldHandler.Christoffel(ig,D1,1,1,2), ManifoldHandler.Christoffel(ig,D1,1,2,2);
                ManifoldHandler.Christoffel(ig,D1,2,1,2), ManifoldHandler.Christoffel(ig,D1,2,2,2)]);
    C = simplify(C);
    DC1 = simplify(cat(4, diff(C, vars(1)), diff(C, vars(2))));
    Schr = ManifoldHandler.RicciScalarChristoffel(ig, C, DC1);
    S = simplify(Schr);    
    
end