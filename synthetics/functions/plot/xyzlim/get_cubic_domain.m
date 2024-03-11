function lim = get_cubic_domain(lat0, lon0, dep0, sidelength)
% compute the x/y/z limits for a cube with its centre at lat0/lon0/dep0

% Conversion between km and degrees
r_earth = 6371.009;
km2lat  = 1/deg2km(1);
lon2km  = r_earth*cos(lat0*pi/180)*pi/180;
km2lon  = 1/lon2km;
%lat2km  = 1/km2lat;

dlat = sidelength/2*km2lat;
dlon = sidelength/2*km2lon;
ddep = sidelength/2;

lim.x = [lat0-dlat lat0+dlat];
lim.y = [lon0-dlon lon0+dlon];
lim.z = [dep0-ddep dep0+ddep];