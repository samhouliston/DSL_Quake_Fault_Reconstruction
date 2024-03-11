function evd = get_eigenval_decomp(cat)

neq = numel(cat.mag);

if neq>3
    
    % Cluster centre
    evd.lat_median = median(cat.lat); 
    evd.lon_median = median(cat.lon); 
    evd.dep_median = median(cat.dep); 
    
    % Cluster geometry
    dlat = cat.lat-evd.lat_median;     % Right sided coordinate system
    dlon = cat.lon-evd.lon_median;
    dz   = cat.dep-evd.dep_median;
    
    dx = deg2km(dlat);
    dy = lon2km(dlon, cat.lat);
    dr = sqrt(dx.^2 +dy.^2 +dz.^2);
    dX = [dx, dy, dz];
    
    % Eigenvalue decomposition
    covmat = cov(dX);
    [V,D]  = eig(covmat);
    
    % Rotate cluster
    dXprime = dX*V;             % Use eigenmatrix as rotation matrix
    dXprime = fliplr(dXprime);  % Matlab returns smallest EV first
    
    % When no. of cat is very small, EVs can be slightly negative
    if min(D(:))<0; evd.Dneg = min(D(:));
                    D=D+1e-12;
    else            evd.Dneg = [];
    end
    
    % Cluster size estimates
    evd.EV = sort(diag(D),'descend')';
    evd.L  = sqrt( evd.EV(1) );
    evd.W  = sqrt( evd.EV(2) );
    evd.T  = sqrt( evd.EV(3) );

    % Save all values for output
    evd.dx      = dx;
    evd.dy      = dy;
    evd.dz      = dz;
    evd.dr      = dr;
    evd.dX      = dX;
    evd.dXprime = dXprime;
    evd.V       = V;
    evd.D       = D;
    

else
    evd.EV         = zeros(3,1);
    evd.lat_median = [];
    evd.lon_median = [];
    evd.dep_median = [];
    evd.dx         = [];
    evd.dy         = [];
    evd.dz         = [];
    evd.dr         = [];
    evd.dX         = [];
    evd.dXprime    = [];
    evd.V          = [];
    evd.D          = [];
    evd.Dneg       = [];
end