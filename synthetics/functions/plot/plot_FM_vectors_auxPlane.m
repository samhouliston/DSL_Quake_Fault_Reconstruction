function plot_FM_vectors_auxPlane(fm, x0, y0, z0)

c = 10;
dx = c/10;

lnw = 3;

% Slip vector
quiver3(x0, y0, z0, ...
    c*fm.s2(1), ...
    c*fm.s2(2), ...
    c*fm.s2(3), ...
    'lineWidth', lnw, ...
    'color','k', ...
    'autoscale', 'off')

text(x0 +c*fm.s2(1) +dx, ...
     y0 +c*fm.s2(2) +dx, ...
     z0 +c*fm.s2(3), ...
     's_2', ...
     'color','k')

% Normal vector    
quiver3(x0, y0, z0, ...
    c*fm.n2(1), ...
    c*fm.n2(2), ...
    c*fm.n2(3), ...
    'lineWidth', lnw, ...
    'color','r', ...
    'autoscale', 'off')

text(x0 +c*fm.n2(1) +dx, ...
     y0 +c*fm.n2(2) +dx, ...
     z0 +c*fm.n2(3), ...
     'n_2', ...
     'color','r')
 

% B-vector    
quiver3(x0, y0, z0, ...
    c*fm.b2(1), ...
    c*fm.b2(2), ...
    c*fm.b2(3), ...
    'lineWidth', lnw, ...
    'color','b', ...
    'autoscale', 'off')

text(x0 +c*fm.b2(1) +dx, ...
     y0 +c*fm.b2(2) +dx, ...
     z0 +c*fm.b2(3), ...
     'b_2', ...
     'color','b')

 
% P-vector    
% (points towards focal origin)
quiver3(...
    x0 +c/2*fm.p2(1), ...
    y0 +c/2*fm.p2(2), ...
    z0 +c/2*fm.p2(3), ...
    -c/2*fm.p2(1), ...
    -c/2*fm.p2(2), ...
    -c/2*fm.p2(3), ...
    'lineWidth', lnw, ...
    'color','m', ...
    'autoscale', 'off')

text(x0 +c/2*fm.p2(1) +dx, ...
     y0 +c/2*fm.p2(2) +dx, ...
     z0 +c/2*fm.p2(3), ...
     'P_2', ...
     'color','m')

% T-vector    
quiver3(x0, y0, z0, ...
    c/2*fm.t2(1), ...
    c/2*fm.t2(2), ...
    c/2*fm.t2(3), ...
    'lineWidth', lnw, ...
    'color','c', ...
    'autoscale', 'off')

text(x0 +c/2*fm.t2(1) +dx, ...
     y0 +c/2*fm.t2(2) +dx, ...
     z0 +c/2*fm.t2(3), ...
     'T_2', ...
     'color','c')