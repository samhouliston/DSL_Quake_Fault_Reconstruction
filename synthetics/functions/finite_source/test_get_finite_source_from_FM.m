clear all

lat0 = 39;
lon0 = 10;
dep0 = 5000;
stk  = 190;
dip  = 60;
mag  = 5;
src  = get_finite_source_from_FM(lat0, lon0, dep0, stk, dip, mag);

figure(1); clf; hold on; grid on; box on;
set(gca,'zDir', 'reverse')
axis ij
xlabel('Lat'); ylabel('Lon'); zlabel('Depth [m]')

% plot3(src.lat5(:), src.lon5(:), src.dep5(:), '-o', ...
%     'markerEdgeColor', [.8 .8 .8], ...
%     'markerFaceColor', [.3 .3 .3])

plot3(lat0, lon0, dep0, 'o', ...
    'markerEdgeColor', [.8 .8 .8], ...
    'markerFaceColor', 'r')

fill3( src.lat5(:), src.lon5(:), src.dep5(:), ones(5,1), ...
    'faceColor', 'b', ...
    'edgeColor', 'k', ...
    'faceAlpha', .3)

nl = 2;
set(gca, 'xlim', [lat0-nl*src.L_deg, lat0+nl*src.L_deg], ...
         'ylim', [lon0-nl*src.L_deg, lon0+nl*src.L_deg], ...
         'zlim', [dep0-nl*src.L    , dep0+nl*src.L])

plot3([-1e5 1e5], [0 0], [0 0],'color', [.7 .7 .7])
plot3([0 0], [-1e5 1e5], [0 0],'color', [.7 .7 .7])
plot3([0 0], [0 0], [-1e5 1e5],'color', [.7 .7 .7])
     
%set(gca,'zDir', 'normal')

%set(gca,'view', [0 90]) % Mapview
set(gca,'view', [-90 90]) % Mapview with 'axis ij'
set(gca,'view', [-100 70]) %
%set(gca,'view', [30 90]) % 
%set(gca,'view', [15 60])
%set(gca,'view', [-20 50])
