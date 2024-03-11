function plot_finite_sources_from_FM_cat(catFM, plt, lim, catH5)
%function plot_finite_sources_from_FM_cat(catFM, lim, sta, F, catH5, finsrc, cmaps)

alpha = .6;

finsrc.plotme = [1 1 1 1 1];
finsrc.plotgrey = true;
%hbm = plot_Amatrice16_basemap(sta, F, catH5, finsrc);
hbm = plot_Amatrice16_basemap(plt);

set(gca,'yDir','reverse', ...
        'zDir','reverse', ...
        'view', [-110 70], ...
        'xlim', lim.x, ...
        'ylim', lim.y, ...
        'zlim', lim.z)

% ax1 = gca;
% colormap(ax1,viridis(20))


% Plot source patches
useme = catFM.lat>=lim.x(1) & catFM.lat<lim.x(2) & ...
        catFM.lon>=lim.y(1) & catFM.lon<lim.y(2) & ...
        catFM.dep>=lim.z(1) & catFM.dep<lim.z(2);

catFM = select_subcat(catFM, useme);

nfm = numel(catFM.mag);

% Colour source patches by rake
crange = [-180 180];
cmap   = plt.cmaps.rb24_2circ;
cnums  = numel(cmap(:,1));
clab   = 'Rake [deg]';

for ifm = 1:nfm
    
    print_iter_nums(ifm,nfm,100)
    
    src = catFM.finsrc{ifm};
    rak = catFM.rak(ifm);
    
    c = sample_colour_from_range(rak, crange, cnums, cmap);

    fill3(src.lat5(:), src.lon5(:), src.dep5(:), ones(5,1), ...
        'faceColor', c, ...
        'edgeColor', 'k', ...
        'faceAlpha', alpha)
end

title(sprintf('Stressdrop = %i MPa', src.stressdrop*1e-6))

cb = colorbar;
cb.Label.String = 'Rake [deg]';
colormap(cmap)
caxis(crange);
cb.Ticks = crange(1):30:crange(2);







%     plot3(src.lat5(:), src.lon5(:), src.dep5(:), '-o', ...
    %         'markerEdgeColor', [.8 .8 .8], ...
    %         'markerFaceColor', [.3 .3 .3])
    
    %     plot3(src.lat0, src.lon0, src.dep0, 'd', ...
    %         'markerEdgeColor', [.8 .8 .8], ...
    %         'markerFaceColor', 'r')
    
% nl = 2;
% set(gca, 'xlim', [src.lat0-nl*src.L_deg, src.lat0+nl*src.L_deg], ...
%          'ylim', [src.lon0-nl*src.L_deg, src.lon0+nl*src.L_deg], ...
%          'zlim', [src.dep0-nl*src.L    , src.dep0+nl*src.L])
% 
% plot3([-1e5 1e5], [0 0], [0 0],'color', [.7 .7 .7])
% plot3([0 0], [-1e5 1e5], [0 0],'color', [.7 .7 .7])
% plot3([0 0], [0 0], [-1e5 1e5],'color', [.7 .7 .7])
     
%set(gca,'zDir', 'normal')
