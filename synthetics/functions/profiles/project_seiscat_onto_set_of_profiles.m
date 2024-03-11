function project_seiscat_onto_set_of_profiles(cat, profiles, rmax, plt)

% Inputs
% cat       hypocentre catalogue
% profiles  set of profiles onto which quakes are projected
% rmax      maximum perpendicular distance from profile
% plt       plotting options and parameters


% Plotting parameters
mec = [.4 .4 .4];
mfc = 'y';
% c   = parula(40); 
% mkc = [.3 .3 .3];
% mfa = .75;
% mea = .6;

nprofiles = numel(profiles.x0);

% Label every ntxt^th profile
ntxt = 10*round( (nprofiles/20)/10 ); 
if ntxt<1; ntxt=1; end
itxt = 0:ntxt:nprofiles;
itxt(1) = 1;
dx    = profiles.dw/10;
dy    = profiles.dw/10;
tstring = sprintf('%s - rmax=%im - nprofiles=%i', ...
    cat.prop.name, ...
        rmax, ...
        nprofiles);

% Prepare movie
if plt.filmme 
    videoFullName = sprintf('%s/new/seiscat_vProfiles_rmax%04dm_spacing%04dm_nprof%04d', ...
        plt.figDir, ...
        rmax, ...
        profiles.dw, ...
        nprofiles);

    V = VideoWriter(videoFullName, 'MPEG-4');
	V.FrameRate = 8; % Frames per second
    V.Quality   = 100;
    open(V);
end









%% Loop over profiles
for iprofile = 1:nprofiles
    
    print_iter_nums(iprofile, nprofiles, 10)
    
    pstring = sprintf('%i / %i profiles of %s - rmax=%im', ...
        iprofile, ...
        nprofiles, ...
        cat.prop.name, ...
        rmax);

    % Compute distance of all points to profile
    [~, r, rS, rD] = project_pointcloud_onto_plane(...
        [cat.n, cat.e, cat.d], ...
        profiles, ...
        iprofile, ...
        0);
    
    % Select hypocentres within rmax of profile
	useme = abs(r)<rmax;
    rS    = rS(useme);
    rD    = rD(useme);
    %n     = sum(useme);
	%cati  = select_subcat(cat,useme);
    
    
    % Plot profile ........................................................
    hf = figure(205); clf; hold on; grid on; box on; axis equal;
    set(hf, 'defaultLegendAutoUpdate','off');
    set(gca,'yDir','reverse')
    xlabel('Along-profile distance [m]')
    ylabel('Depth [m]')
    title(pstring)
    
    % Plot along strike distance vs absolute depth
    plot(rS, rD+profiles.z0(iprofile), ...
         '.', ...
         'color', 'k', ...
         'DisplayName', 'Hypocentres');

    % Plotting parameters
    L    = profiles.length;
    Dtop = profiles.z0(iprofile);
    Dbot = profiles.z0(iprofile) + profiles.depth;
    
    set(gca, 'xlim', [0 L], ...
             'ylim', [Dtop Dbot])
	
    
    % Add overview map inlet ..............................................
    ha = axes; box on; hold on;
    axis equal
    ha.Position(1) = ha.Position(1) +0.05;
    ha.Position(2) = ha.Position(2) +0.05;
    ha.Position(3) = 0.2;
    ha.Position(4) = 0.2; %ha.Position = [.14 .3 .12 .15];
    %ha.XTickLabel = '';
    %ha.YTickLabel = '';
    plot(cat.n, cat.e, '.k')
    xlabel('')
    ylabel('')
    
    %set(gca,'view', plt.view_angle, ...
    set(gca,'view', [-90 90], ...
            'yDir','reverse', ...
            'zDir','reverse')

	% Zoom out to see all profiles in map
    set_bounding_box([ [profiles.x0; profiles.xE], ...
                       [profiles.y0; profiles.yE] ]);

    % Show all profile lines
    plot([profiles.x0 profiles.xE]', ...
         [profiles.y0 profiles.yE]', ...
        '-y')

    % Show current, single profile line as arrow
    quiver(profiles.x0(iprofile), ...
           profiles.y0(iprofile), ...
           profiles.xE(iprofile) -profiles.x0(iprofile), ...
           profiles.yE(iprofile) -profiles.y0(iprofile), ...
           'lineWidth', 1, ...
            'color', 'r', ...
            'autoscale', 'off');
    
    text(profiles.x0(iprofile) -dx, ...
         profiles.y0(iprofile) -dy, ...
         sprintf('%i', iprofile), ...
         'HorizontalAlignment', 'right');

    
    
    
    % hl = legend;
    % hl.Location = 'southWest';
    % hl.TextColor = [.2 .2 .2];
    
    
    
    
    

    
    % Save frame in movie .................................................
    if plt.filmme
        frame = getframe(gcf);
        writeVideo(V, frame);
    end
    
    % Print figure ........................................................
    if plt.printme
        set(gcf,'paperposition',[0 0 16 12])
        
        figName = sprintf('%s/new/FM_profiles_%04d', ...
            plt.figDir, ...
            iprofile);
        
        print('-dpng', figName)
    end
end

if plt.filmme; close(V);
end






































%% APPENDIX


% %% Plot profile overview map 
% nreps = 5;
% for irep = 1:nreps
%     
%     figure(7322); clf;
%     hbm = plot_Amatrice16_basemap(plt);
%     axis equal
%     set(gca,'view', [-90 90], ...
%         'xlim', profiles.xlm, ...
%         'ylim', profiles.ylm, ...
%         'zlim', profiles.zlm)
% 
%     plot3(profiles.x0, profiles.y0, profiles.z0, 's', ...
%         'markerFaceColor', mfc, ...
%         'markerEdgeColor', mec, ...
%         'markerSize',5)
%     
%     plot3(profiles.xE, profiles.yE, profiles.zE, 's', ...
%         'markerFaceColor', mfc, ...
%         'markerEdgeColor', mec, ...
%         'markerSize',5)
%     
%     for iprofile = 1:nprofiles
%         
%         if ismember(iprofile, itxt)
%             
%             plot([profiles.x0(iprofile) profiles.xE(iprofile)], ...
%                  [profiles.y0(iprofile) profiles.yE(iprofile)], ...
%                   '-k')
%             
%             text(profiles.x0(iprofile) -dlat, ...
%                  profiles.y0(iprofile) -dlon, ...
%                  sprintf('%i', iprofile), ...
%                  'HorizontalAlignment', 'right');
%         end
%     end
%     
%     l_cross = round(profiles.length/1e3/5);
%     plot_dlat_dlon_cross(l_cross)
%     
%     
%     
%     if plt.filmme
%         frame = getframe(gcf);
%         writeVideo(V, frame);
%     end
% end

% scatter3(cati.x, cati.y, cati.z, fl*cati.length, ...
%     r(useme), ...
%     'filled', ...
%     'markerEdgeColor','k', ...
%     'markerEdgeAlpha', mea, ...
%     'markerFaceAlpha', mfa);

% %scatter3(catFMi.x, catFMi.y, catFMi.z, fl*catFMi.length, ...
% scatter3(catFMi.x, catFMi.y, -.1*ones(size(catFMi.z)), fl*catFMi.length, ...
%     r_FM(useme_FM), ...
%     'filled', ...
%     'markerEdgeColor','r', ...
%     'lineWidth', 2, ...
%     'markerEdgeAlpha', mea, ...
%     'markerFaceAlpha', mfa);



    %     set_bounding_box([[profiles.x0; profiles.xE], ...
    %                       [profiles.y0; profiles.yE]], ...
    %                      100);
    %     set(gca,'zlim', [-1 20])
    

%     [r, rS, rD] = get_projection_of_points_onto_known_profile(points, ...
%                       profiles, ...
%                       iprofile, ...
%                       0);
                  
    


    %     mksize = 30;
    %     mfa    = .75;
    %     mea    = .6;
    %     fl     = 1;
    % 
    %     cmaps = plt.cmaps;


    %     scatter(rS(useme), rD(useme), fl*cati.length, ...
    %         r(useme), ...
    %         'filled', ...
    %         'markerEdgeColor','k', ...
    %         'markerEdgeAlpha', mea, ...
    %         'markerFaceAlpha', mfa);
    %     hc = colorbar;
    %     colormap(hs1, cmaps.kagan_angle)
    %     caxis([-rmax rmax])
    %     hc.Label.String = 'Fault-perpendicular distance [m]';

    

%     % 2. Seismicity 3D map ................................................
%     subplot(hs2);
%     hbm = plot_Amatrice16_basemap(plt);
%     
%     plot3(P0(1),  P0(2),  P0(3),  'xr', 'markerSize',20, 'lineWidth',3)
%     
%     set_bounding_box([cat.x, cat.y, cat.z], 99.8);
%     
%     set(gca,'zlim', [zmin zmax], ...
%             'view', [-90 90], ... %'view', [-70 65], ...
%             'yDir','reverse', ...
%             'zDir','reverse')