function setup_Carlsson_KB(results_dir)

    DCT = load(sprintf('%s/DCTParams.mat', results_dir));

    fprintf('Generating H\n');

    [Hx, Hy] = meshgrid(-1:1, 1:-1:-1);
    Hx = Hx(:);
    Hy = Hy(:);

    fprintf('Generating K\n');

    syms theta phi
    a = cos(theta);
    b = sin(theta);
    c = cos(phi);
    d = sin(phi);

    P = @(x,y) c*(a*x+b*y)^2 + d*(a*x+b*y);

    P_of_H = simplify(arrayfun(@(x,y) P(x,y), Hx, Hy));   
    Klein_R9_fcn = matlabFunction(P_of_H', 'Vars', [theta, phi]);    
    Dnorm = simplify(sqrt(sum((DCT.D*P_of_H).*P_of_H,1)));    
    P_of_H_norm = simplify(DCT.lambda*DCT.A'*(P_of_H./Dnorm));
    Klein_S7_fcn = matlabFunction(P_of_H_norm', 'Vars', [theta, phi]);
    clear DCT;

    fprintf('Generating S Expression\n');

    g = ManifoldHandler.get_metric(P_of_H_norm, theta, phi); 
    det_g = simplify(det(g));
    det_g_fcn = matlabFunction(det_g, 'Vars', [theta, phi]);

    Schr = ManifoldHandler.get_scal_curv(g, theta, phi);
    Schr_fcn = matlabFunction(Schr, 'Vars', [theta, phi]);
    
    save(sprintf('%s/CarlssonKB.mat', results_dir));
    
end
