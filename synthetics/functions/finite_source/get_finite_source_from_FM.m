function src = get_finite_source_from_FM(lat0, lon0, dep0, strike, dip, mag, stressdrop)
% Compute square-shaped finite source centered around hypocentre.

% Conversion between km and degrees
r_earth = 6371.009;
km2lat  = 1/deg2km(1);
lon2km  = r_earth*cos(lat0*pi/180)*pi/180;
km2lon  = 1/lon2km;
%lat2km  = 1/km2lat;


% Rectangular source model used in Meier et al., 2014, JGR
% estimate source patch dimension L of a square source with the relation
% for strike-slip faults [Knopoff, 1958; Scholz, 2002]
% where L is the fault dimension in meters.
% estimate average slip u over the fault area L , use def of seismic moment
% [Aki and Richards, 2002]

%stressdrop = 1e6;  % [MPa]
shearmod   = 30e9; % [GPa]

M0 = magnitude2moment(mag);
L  = ( 2*M0./(pi*stressdrop) ).^(1/3); % [m]
D  = M0./(shearmod*L.^2);     % [m]

% Change to radians
dipdir = (strike-90)*pi/180;
dip    = dip        *pi/180;
strike = strike     *pi/180;

% Compute 4 corners of square with centre at [x0, y0, z0]
rs = [-L/2; +L/2; +L/2; -L/2];
rd = [-L/2; -L/2; +L/2; +L/2];

dx =  rs*cos(strike) +rd*cos(dip)*cos(dipdir);
dy =  rs*sin(strike) +rd*cos(dip)*sin(dipdir);
dz = -rd*sin(dip);

dx_lat = dx/1000*km2lat;
dy_lon = dy/1000*km2lon;

src.lat4 = lat0 + dx_lat;
src.lon4 = lon0 + dy_lon;
src.dep4 = dep0 + dz/1000;

% Make closed surface, for plotting
src.lat5 = [src.lat4; src.lat4(1)];
src.lon5 = [src.lon4; src.lon4(1)];
src.dep5 = [src.dep4; src.dep4(1)];

src.lat0       = lat0;
src.lon0       = lon0;
src.dep0       = dep0;
src.stk        = strike*180/pi;
src.dip        = dip*180/pi;
src.dipdir     = (src.stk-90);
src.L          = L;
src.L_deg      = L/1000*km2lon;
src.D          = D;
src.stressdrop = stressdrop;
src.shearmod   = shearmod;

% [RS,RD] = meshgrid(rs,rd);
% x = x0 + RS*cos(strike) + RD*cos(dip)*cos(dipdir);
% y = y0 + RS*sin(strike) + RD*cos(dip)*sin(dipdir);
% z = z0 + RD*sin(dip);