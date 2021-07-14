classdef (Sealed=true) ManifoldHandler
        
    methods (Static)
        
        Cijk = Christoffel(ig, D1, i, j, k);
        S = RicciScalarChristoffel(ig, C, DC1);
        
        [X, vars] = uniform_sample(X, N_points, bounds)
        [X, vars] = rejection_sample(X, N_points, target_density, bounds, skew)
        g = get_metric(X);
        S = get_scal_curv(X);
        
        %% Hypersphere - Sn
        
        function X = hypersphere_S2(R)
            syms theta phi
            X = [R.*cos(theta).*cos(phi); R.*cos(theta).*sin(phi); R.*sin(theta)];
        end
        
        function X = sample_hypersphere(N, dim, R)            
            if (nargin < 3)
                R = 1;
            end
            % Equivalent to rejection sampling from the metric    
            X = randn(N, dim+1);
            X = R.*X./sqrt(sum(X.^2, 2));
        end
        
        function S = hypersphere_curvature(dim, R)
            if (nargin < 2)
                R = 1;
            end
            S = dim.*(dim-1)./R;
        end
        
        %% One-Sheet Hyperboloid - H2
        
        function X = one_sheet_hyperboloid(a, b, c)
            syms theta u
            X = [a*sqrt(u.^2+1).*cos(theta); b*sqrt(u.^2+1).*sin(theta); c*u];
        end
        
        function [X, coords] = sample_one_sheet_hyperboloid(N, a, b, c)
            X = ManifoldHandler.one_sheet_hyperboloid(a, b, c);
            [X, coords] = ManifoldHandler.uniform_sample(X, N, [0, 2*pi; -2 2]);
        end
          
        function S = one_sheet_hyperboloid_curvature(z, a, b, c)          
            S = ManifoldHandler.get_scal_curv(ManifoldHandler.one_sheet_hyperboloid(a, b, c));
            S = matlabFunction(S, 'Vars', symvar(S));
            S = arrayfun(@(z) S(z/c), z);
            % S = -2*c^2./((a^2+c^2)*(z/c).^2 + c^2).^2;
        end
        
        %% Torus - T2
            
         function X = torus(R, r)
            syms inner outer
            X = [(R+r.*cos(inner)).*cos(outer); (R+r.*cos(inner)).*sin(outer); r.*sin(inner)];
        end
            
        function [X, angles] = sample_torus(N, R, r)
            X = ManifoldHandler.torus(R, r);
            [X, angles] = ManifoldHandler.uniform_sample(X, N, [0, 2*pi; 0, 2*pi]);
        end
        
        function S = torus_curvature(R, r, theta)
            assert(R > r);
            S = ManifoldHandler.get_scal_curv(ManifoldHandler.torus(R, r));
            S = matlabFunction(S, 'Vars', symvar(S));
            S = arrayfun(@(t) S(t), theta);
            % S = 2*cos(inner)./(r*(R+r*cos(inner)));
        end
        
        function X = augment_dimension(X, dim, method)
            [N, d] = size(X);
            if (dim == d)
                return;
            end
            assert(dim > d);
            if (strcmp(method, 'zeropad'))
                X = [X zeros(N, dim-d)];
            elseif (strcmp(method, 'pca'))
                nullspace = randn(dim, dim-d);
                basis = null(nullspace');                
                X = X*basis';
            end
        end
        
    end
end