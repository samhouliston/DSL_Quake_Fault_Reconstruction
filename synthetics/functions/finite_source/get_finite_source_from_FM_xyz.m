function src = get_finite_source_from_FM_xyz(x0, y0, z0, strike, dip, mag)
%function get_finite_source_from_FM(lat0, lon0, dep0, strike, dip, mag)





%% Model used in Meier et al., 2014, JGR
% estimate source patch dimension L of a square source with the relation
% for strike-slip faults [Knopoff, 1958; Scholz, 2002]
% where L is the fault dimension in meters.
% estimate average slip u over the fault area L , use def of seismic moment
% [Aki and Richards, 2002]
stressdrop = 3e6;  % [MPa]
shearmod   = 30e9; % [GPa]

M0 = magnitude2moment(mag);
L  = ( 2*M0./(pi*stressdrop) ).^(1/3); % [m]
D  = M0./(shearmod*L.^2);     % [m]

% rs      = -L/2 : L : L/2;  % Distance along azi direction [m]
% rd      = -L/2 : L : L/2;  % Distance along dip direction [m]

% Change to radians
dipdir = (strike-90)*pi/180;
dip    = dip        *pi/180;
strike = strike     *pi/180;

% Compute 4 corners of square with centre at [x0, y0, z0]
rs = [-L/2; ...
      +L/2; ...
      +L/2; ...
      -L/2];

rd = [-L/2; ...
      -L/2; ...
      +L/2; ...
      +L/2];

src.x4 = x0 + rs*cos(strike) + rd*cos(dip)*cos(dipdir);
src.y4 = y0 + rs*sin(strike) + rd*cos(dip)*sin(dipdir);
src.z4 = z0 + rd*sin(dip);

% Make closed surface
src.x5 = [src.x4; src.x4(1)];
src.y5 = [src.y4; src.y4(1)];
src.z5 = [src.z4; src.z4(1)];

src.stk = strike;
src.dip = dip;
src.L   = L;
src.D   = D;
src.stressdrop = stressdrop;
src.shearmod   = shearmod;

% [RS,RD] = meshgrid(rs,rd);
% x = x0 + RS*cos(strike) + RD*cos(dip)*cos(dipdir);
% y = y0 + RS*sin(strike) + RD*cos(dip)*sin(dipdir);
% z = z0 + RD*sin(dip);


% % Conversion between km and degrees
% r_earth = 6371.009;
% km2lat  = 1/deg2km(1);
% lat2km  = 1/km2lat;
% lon2km  = r_earth*cos(lat0*pi/180)*pi/180;
% km2lon  = 1/lon2km;




% % High resolutioon fault grid
% 
% 
% % Source geometry
% dx     = L/2*cos(dip);        % [km]
% dz     = L/2*sin(dip);        % [km]
% d1_lon = dx*cos(strike);      % [km]
% d1_lat = dx*sin(strike);      % [km]
% d2_lon = L/2*sin(strike);     % [km]
% d2_lat = L/2*cos(strike);     % [km]
% 
% d1_lon_deg = d1_lon*km2lon;   % [deg]
% d1_lat_deg = d1_lat*km2lat;   % [deg]
% d2_lon_deg = d2_lon*km2lon;   % [deg]
% d2_lat_deg = d2_lat*km2lat;   % [deg]
% 
% slat = lat0 +d1_lat_deg -d2_lat_deg; % [deg]
% slon = lon0 -d1_lon_deg -d2_lon_deg; % [deg]
% sdep = dep0 -dz;                     % [km]
% 
% % if (sdep<0)                        % code cannot handle events
% %     sdep = 0;                      % with depths < 0
% % end
% 
% % --> For latitude, longitude and depth, save one column per source as:
% % Slat =   [upper start corner          src2      src3     ...
% %           upper end corner            ...       ...      ...
% %           lower start corner          ...       ...      ...
