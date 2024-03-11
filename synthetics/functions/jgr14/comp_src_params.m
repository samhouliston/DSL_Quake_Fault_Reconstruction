function [lat_ref,lon_ref,z_ref,u_strike,u_dip] = comp_src_params(lat,lon,depth,strike,dip,rake,l_rect,u_rect)
% computes coordinates of 1 upper edge of a square fault centered around a
% hypocenter at lat/lon

% UNITS
% input:  angles [RAD], u_rect [m], l_rect & depth [km]
% output: z_ref [km], u_strike & u_dip [m]

% Convert km to degrees
km2lat = 1/111.2;
[one_deg_lon_in_km] = lon2km_smallCircle(1,lat);
km2lon = 1/one_deg_lon_in_km;

d_hor       = l_rect/2*cos(dip);        % [km]
dz          = l_rect/2*sin(dip);        % [km]

d1_lon      = d_hor*cos(strike);        % [km]
d1_lon_deg  = d1_lon*km2lon;            % [deg]
d2_lon      = l_rect/2*sin(strike);     % [km]
d2_lon_deg  = d2_lon*km2lon;            % [deg]

d1_lat      = d_hor*sin(strike);        % [km]
d1_lat_deg  = d1_lat*km2lat;            % [deg]
d2_lat      = l_rect/2*cos(strike);     % [km]
d2_lat_deg  = d2_lat*km2lat;            % [deg]

lat_ref = lat+d1_lat_deg-d2_lat_deg;    % [deg]
lon_ref = lon-d1_lon_deg-d2_lon_deg;    % [deg]
z_ref   = depth-dz;                     % [km]

    if (z_ref<0)                        % code cannot handle events
        z_ref = 0;                      % with depths < 0
    end

u_strike      = u_rect*cos(rake);       % [m]
u_dip         = -u_rect*sin(rake);      % [m]