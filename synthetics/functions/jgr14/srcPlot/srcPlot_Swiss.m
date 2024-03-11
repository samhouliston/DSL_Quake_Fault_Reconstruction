function [Slat,Slon,Sz] = srcPlot_Swiss(lat,lon,depth,strike,dip,src_length)
% Computes 3D coordinates of the four corners of a square fault centered 
% around a hypocenter at lat/lon/depth
% UNITS: input angles in [deg], src_length [km], depth [km]

Neq  = length(lat);
Slat = zeros(4,Neq);
Slon = zeros(4,Neq);
Sz   = zeros(2,Neq);

figure
title('Square fault patches')
hold on
set(gca, 'ZDir', 'reverse')
grid on
xlabel('X','FontSize',12,'FontWeight','b')
ylabel('Y','FontSize',12,'FontWeight','b')
zlabel('Z','FontSize',12,'FontWeight','b')
%axis([min(lon)-0.2 max(lon)+0.2 min(lat)-0.2 max(lat)+0.2 0 20])

az = -30;             % Perspective
el = 10;
view(az, el);

for i=1:1:Neq
    
    % deg2rad - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    strike(i) = strike(i)*pi/180;
    dip(i)    = dip(i)*pi/180;
    % -------------------------------------------------------------------------
    
    d_hor       = src_length(i)/2*cos(dip(i));        % [km]
    dz          = src_length(i)/2*sin(dip(i));        % [km]
    
    d1_lon      = d_hor*cos(strike(i));               % [km]
    d2_lon      = src_length(i)/2*sin(strike(i));     % [km]
    
    d1_lat      = d_hor*sin(strike(i));               % [km]
    d2_lat      = src_length(i)/2*cos(strike(i));     % [km]
    
    slat   = lat(i)+d1_lat-d2_lat;    % [deg]
    slon   = lon(i)-d1_lon-d2_lon;    % [deg]
    sz     = depth(i)-dz;                     % [km]
    
    if (sz<0)                        % code cannot handle events
        sz = 0;                      % with depths < 0
    end
        
    % --> For latitude, longitude and depth, save one column per source as:
    % Slat =   [upper start corner          src2      src3     ...
    %           upper end corner            ...       ...      ...
    %           lower start corner          ...       ...      ...
    %           lower end corner            ...       ...      ... ]
    
    % Compute upper end corner
    elat   = slat + src_length(i)*cos(strike(i));
    elon   = slon + src_length(i)*sin(strike(i));
    
    % Upper start- and end-corner
    Slat(1,i) = slat;             % Write pairs of point into matrices of
    Slat(2,i) = elat;             % dimensions 2 x Neq to combine coordinates
    Slon(1,i) = slon;             % with lines pairwise
    Slon(2,i) = elon;
    Sz(1,i)   = sz;
    
    % Compute lower boundary coordinates of fault patch
    ds1 = lat(i) - slat;
    ds2 = lon(i) - slon;
    ds3 = lat(i) - elat;
    ds4 = lon(i) - elon;
    
    slat_low = elat + 2*ds3;
    slon_low = elon + 2*ds4;
    elat_low = slat + 2*ds1;
    elon_low = slon + 2*ds2;
    
    sz_low = depth(i) + dz;                     % [km]

    
    % Lower start- and end-corner
    Slat(3,i) = slat_low;             % Write pairs of point into matrices of
    Slat(4,i) = elat_low;             % dimensions 2 x Neq to combine coordinates
    Slon(3,i) = slon_low;             % with lines pairwise
    Slon(4,i) = elon_low;
    Sz(2,i)   = sz_low;
    
    %% Plot
    %gcf
    
    %plot3(lon(i),lat(i),depth(i),'*r')   % Hypocenter
    %plot3(slon,slat,sz,'*r')             % Reference point
    plot3([Slon(1,i) Slon(2,i)],[Slat(1,i) Slat(2,i)],[Sz(1,i) Sz(1,i)])
    plot3([Slon(3,i) Slon(4,i)],[Slat(3,i) Slat(4,i)],[Sz(2,i) Sz(2,i)])
    plot3([Slon(1,i) Slon(3,i)],[Slat(1,i) Slat(3,i)],[Sz(1,i) Sz(2,i)])
    plot3([Slon(2,i) Slon(4,i)],[Slat(2,i) Slat(4,i)],[Sz(1,i) Sz(2,i)])
    
    plot3([Slon(1,i) Slon(4,i)],[Slat(1,i) Slat(4,i)],[Sz(1,i) Sz(2,i)],'k')
    plot3([Slon(2,i) Slon(3,i)],[Slat(2,i) Slat(3,i)],[Sz(1,i) Sz(2,i)],'k')
    
%     x = linspace(slon,elon,3);
%     y = linspace(slat,elat,3);
%     z = [sz, sz, sz; sz+(sz_low-sz)/2, sz+(sz_low-sz)/2, sz+(sz_low-sz)/2; sz_low,sz_low,sz_low]
%     
%     figure(2)
%     mesh(x,y,z)
%     set(gca, 'ZDir', 'reverse')
title(['event nr. ',num2str(i)])
    
end
