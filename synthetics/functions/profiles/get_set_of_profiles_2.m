function p = get_set_of_profiles_2(p)

% Same as get_set_of_profiles.m but instead of specifying strike and dip,
% one specifies the n/e/d coordinates of the first and last profile
% starting point. 


p.x0 = linspace(p.x1, p.x2, p.nprof);
p.y0 = linspace(p.y1, p.y2, p.nprof);
p.z0 = linspace(p.z1, p.z2, p.nprof);

% % Starting points of all profiles
% ww   = (0 : p.dw : p.wmax)';
% p.x0 = p.x0 + ww*cosd(p.strike -90);
% p.y0 = p.y0 + ww*sind(p.strike -90);
% p.z0 = repmat(p.z0, size(p.x0) );
% 
% % End points of all profiles
% p.xE = p.x0 +p.length*cosd(p.strike);
% p.yE = p.y0 +p.length*sind(p.strike);
% p.zE = p.z0;



%% Normal vector

% Calculate normal vector of the plane with strike and dip
nx = -sind(p.dip) * sind(p.strike);
ny =  sind(p.dip) * cosd(p.strike);
nz = -cosd(p.dip);
p.normal_vector = [nx, ny, nz]; % normal vector

% Find d so that plane goes through point P0
p.d  = -nx*p.x0 -ny*p.y0 -nz*p.z0; %d  = -nx*P0(1) -ny*P0(2) -nz*P0(3);

% Xnx = -sind(dip) * sind(strike); % Normal vector definition for Stein &
% Xny = -sind(dip) * cosd(strike); % Wysession does not work, presumably  
% Xnz =  cosd(dip);                % because we use a North/East/Down 
                                   % coordinate system


% % Compute suitable map boundaries 
% p.xlm = [min([p.x0; p.xE])-p.dw; ...
%          max([p.x0; p.xE])+p.dw];
% p.ylm = [min([p.y0; p.yE])-p.dw; ...
%          max([p.y0; p.yE])+p.dw];        
% p.zlm = [p.z0(1) p.z0(1)+1.1*p.depth];
    





%% Visualise
if p.plotme
    
    mec = [.2 .2 .2];
    %mea = .8;
    
    if p.newfig
        figure(7321); clf; hold on; grid on; box on;
        axis equal;
        set(gca,'view', [-90 90])
        set(gca,'xlim', p.xlm, ...
            'ylim', p.ylm, ...
            'zlim', p.zlm)
        %set_bounding_box([cat.x, cat.y, cat.z], 99.5)
    end
    
    nprofiles = numel(p.x0);
    tstring = sprintf('%i profiles, every %im, with strike=%i°', ...
        nprofiles, ...
        p.dw, ...
        p.strike);
    title(tstring)
    
    % Start points
    plot3(p.x0, p.y0, p.z0, 'o', ...
        'markerFaceColor', 'y', ...
        'markerEdgeColor', mec, ...
        'lineWidth', .1, ...
        'markerSize',8)
    
    % End points
%     plot3(p.xE, p.yE, p.zE, 'o', ...
%         'markerFaceColor', 'y', ...
%         'markerEdgeColor', mec, ...
%         'lineWidth', .1, ...
%         'markerSize',8)
%     
%     % All profile lines
%     plot3([p.x0 p.xE]', ...
%          [p.y0 p.yE]', ...
%          [p.z0 p.z0]', ...
%          '-', ...
%          'color', [.4 .4 .4])
% 
%     % Single profile as arrow
%     iprofile = round(numel(ww)/2);
%     quiver3(p.x0(iprofile), ...
%             p.y0(iprofile), ...
%             p.z0(iprofile), ...
%             p.xE(iprofile) -p.x0(iprofile), ...
%             p.yE(iprofile) -p.y0(iprofile), ...
%             p.zE(iprofile) -p.z0(iprofile), ...
%            'lineWidth', 2, ...
%             'color', 'r', ...
%             'autoscale', 'off');
%     
    %     plot3([p.x0(iprofile) p.xE(iprofile)], ...
    %         [p.y0(iprofile) p.yE(iprofile)], ...
    %         [p.z0(iprofile) p.zE(iprofile)], ...
    %         '-', 'color', mec)
end

% % Zoom out to see all profiles in map
% set_bounding_box([ [p.x0; p.xE], ...
%                    [p.y0; p.yE] ]);