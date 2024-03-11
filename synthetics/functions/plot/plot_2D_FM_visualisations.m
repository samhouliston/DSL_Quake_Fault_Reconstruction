function hf = plot_2D_FM_visualisations(fm)

hf = figure(202); clf; whitebg([.5 .5 .5])
s1 = subplot(1,3,1); hold on;
s2 = subplot(1,3,2); hold on; 
s3 = subplot(1,3,3); hold on; 

% Beachball
subplot(s1);
beachball(fm.stk, fm.dip, fm.rak)

% Stereonet
subplot(s2);
stereonet
[azi1, dip1] = norm_ned2azdip(fm.n1); % Plot both normal vectors
[azi2, dip2] = norm_ned2azdip(fm.n2);
[x_s1, y_s1] = stereoplot_point_cords(azi1, dip1);
[x_s2, y_s2] = stereoplot_point_cords(azi2, dip2);
plot(x_s1, y_s1, 'o', 'markerFaceColor', 'r', 'markerEdgeColor', [.4 .4 .4])
plot(x_s2, y_s2, 'o', 'markerFaceColor', 'b', 'markerEdgeColor', [.4 .4 .4])
text(x_s1+.1, y_s1, 'n_1', 'color','r')
text(x_s2+.1, y_s2, 'n_2', 'color','b')

% Hudson Net
% (I think this may assume a ENU coordinate system, instead of the NED used here)
subplot(s3);
M = fm2moment_tensor(fm.s1, fm.n1); % Compute moment tensor equivalent to FM
drawhudsonnet(M);
title('Hudson Net')