function evd = get_eigenval_decomp_for_reference_events(catH5, evd)


% Select reference events, i.e. events outside cluster but still nearby
dkm    = 3*evd.L;      % .. within 2x largest cluster dimension
ddeg   = km2deg(dkm); 
useme  = catH5.lat >= evd.lat_median -ddeg & ...
         catH5.lat <  evd.lat_median +ddeg & ...
         catH5.lon >= evd.lon_median -ddeg & ...
         catH5.lon <  evd.lon_median +ddeg & ...
         catH5.dep >= evd.dep_median -dkm  & ...
         catH5.dep <  evd.dep_median +dkm;

catH5r = select_subcat(catH5, useme);

% Reference cluster geometry
dlat = catH5r.lat-evd.lat_median;     % Right sided coordinate system
dlon = catH5r.lon-evd.lon_median;
dz   = catH5r.dep-evd.dep_median;

dx = deg2km(dlat);
dy = lon2km(dlon, catH5r.lat);
dX = [dx, dy, dz];

% Rotate cluster, using saved EV decomposition from clustered events
dXprime = dX*evd.V;             % Use eigenmatrix as rotation matrix
dXprime = fliplr(dXprime);      % Matlab returns smallest EV first

evd.dX_ref      = dX;
evd.dXprime_ref = dXprime;

evd.iH5_ref = find(useme);