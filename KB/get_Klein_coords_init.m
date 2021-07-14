function [theta, phi] = get_Klein_coords_init(v, C, alpha)

    assert(size(v,2)==8);    
    assert(length(alpha)==C.N.total);
    
    % tan(theta) = v1/v2 = (dby)/(dax) = ((sin phi)*(sin theta))/((sin phi)*(cos theta))    
    
    vp = zeros(size(v,1),4);
    vp(:,1) = (v(:,1)-sqrt(2)/2*v(:,6))/alpha(C.init_alpha_ind.y);
    vp(:,2) = (v(:,2)-sqrt(2)/2*v(:,7))/alpha(C.init_alpha_ind.x);
    vp(:,3) = (v(:,3)-v(:,8))          /alpha(C.init_alpha_ind.y2);
    vp(:,4) = (v(:,4)-v(:,8))          /alpha(C.init_alpha_ind.x2);
    
    theta = atan(vp(:,1)./-vp(:,2));   % ambiguity modulo pi    
    theta(theta<0) = theta(theta<0)+pi; % assume theta = [0,pi] so no loss
    theta(isnan(theta)) = pi/2;
   
    % tan(phi) = sqrt(v1^2+v2^2)/(v3 + v4) = sin phi / cos phi 
    phi = atan(sqrt(vp(:,1).^2+vp(:,2).^2)./(vp(:,3)+vp(:,4))); % ambiguity modulo pi    
    phi(phi<0) = phi(phi<0)+pi;
    phi(isnan(phi)) = pi/2;
    
    [~, idx] = sort(abs(v(:,[1, 2, 6 7])),2,'descend');
    idx = idx(:,1);
    swap = false(size(idx));
    
    % only care about signs so dropping (.)^2 terms
    
    fcns = {@(t,p) sin(p).*sin(t); 
            @(t,p) -sin(p).*cos(t);
            @(t,p) cos(p);
            @(t,p) cos(p);
            @(t,p) -cos(p).*cos(t).*sin(t);
            @(t,p) sin(p).*sin(t);
            @(t,p) -sin(p).*cos(t);
            @(t,p) cos(p)};        
            
    fcns = fcns([1 2 6 7]);
    
    for b = 1:4
        mask = find(idx==b);
        swap(mask(sign(v(mask,b)) ~= sign(fcns{b}(theta(mask), phi(mask))))) = true;
    end
        
    phi(swap) = 2*pi-phi(swap);
    
end