function evd = get_eigenval_decomp_NED(cat)
% Estimate cluster geometry via a eigenvalue decomposition / principal 
% component analyis, in a North/East/Down coordinate system
% Men-Andrin Meier, 26/5/2023

neq = numel(cat.n);

if neq>3
    
    % Cluster centre
    evd.n_median = median(cat.n); 
    evd.e_median = median(cat.e); 
    evd.d_median = median(cat.d); 
    
    % Cluster geometry
    dn = cat.n-evd.n_median;     % Right sided coordinate system
    de = cat.e-evd.e_median;
    dd = cat.d-evd.d_median;
    
    dr = sqrt(dn.^2 +de.^2 +dd.^2);
    dX = [dn, de, dd];
    
    % Eigenvalue decomposition
    covmat = cov(dX);
    [V,D]  = eig(covmat);
    % D: diagonal matrix of eigenvalues
    % V: full matrix whose columns are the corresponding eigenvectors  
    % so that A*V = V*D.
    
    % Matlab returns smallest EV first; change to largest one first
    V = fliplr(V); 
    
    % Rotate cluster; use eigenmatrix as rotation matrix
    dXprime = dX*V;   
    
    % Pre 26/5/2023:
    %     % Rotate cluster
    %     dXprime = dX*V;             % Use eigenmatrix as rotation matrix
    %     dXprime = fliplr(dXprime);  % Matlab returns smallest EV first
    
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
    evd.dn      = dn;
    evd.de      = de;
    evd.dd      = dd;
    evd.dr      = dr;
    evd.dX      = dX;
    evd.dXprime = dXprime;
    evd.V       = V;
    evd.D       = D;
    

else
    evd.EV         = zeros(3,1);
    evd.n_median = [];
    evd.e_median = [];
    evd.d_median = [];
    evd.dn         = [];
    evd.de         = [];
    evd.dd         = [];
    evd.dr         = [];
    evd.dX         = [];
    evd.dXprime    = [];
    evd.V          = [];
    evd.D          = [];
    evd.Dneg       = [];
end