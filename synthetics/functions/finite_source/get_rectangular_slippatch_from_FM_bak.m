function src = get_rectangular_slippatch_from_FM(...
    x0, y0, z0, ...
    strike, dip, ...
    mag, ...
    stressdrop)
% Compute square-shaped finite source centered around hypocentre.

% Rectangular source model used in Meier et al., 2014, JGR
% estimate source patch dimension L of a square source with the relation
% for strike-slip faults [Knopoff, 1958; Scholz, 2002]
% where L is the fault dimension in meters.
% estimate average slip u over the fault area L , use def of seismic moment
% [Aki and Richards, 2002]

% Use the following commands and script for plotting:
% plot_rectangular_slippatch(src, x0, y0, z0, colour)

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

dx = rs*cos(strike) +rd*cos(dip)*cos(dipdir);
dy = rs*sin(strike) +rd*cos(dip)*sin(dipdir);
dz =                -rd*sin(dip);

src.x4 = x0 +dx;
src.y4 = y0 +dy;
src.z4 = z0 +dz;

% Make closed surface, for plotting with fill3.m
src.x5 = [src.x4; src.x4(1)];
src.y5 = [src.y4; src.y4(1)];
src.z5 = [src.z4; src.z4(1)];

% Store everything in output structure 
src.x0         = x0;
src.y0         = y0;
src.z0         = z0;
src.mag        = mag;
src.stk        = strike*180/pi;
src.dip        = dip*180/pi;
src.dipdir     = (src.stk-90);
src.length     = L;
src.slip       = D;
src.stressdrop = stressdrop;
src.shearmod   = shearmod;

% [RS,RD] = meshgrid(rs,rd);
% x = x0 + RS*cos(strike) + RD*cos(dip)*cos(dipdir);
% y = y0 + RS*sin(strike) + RD*cos(dip)*sin(dipdir);
% z = z0 + RD*sin(dip);