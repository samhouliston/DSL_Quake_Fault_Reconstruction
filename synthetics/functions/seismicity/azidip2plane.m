function [x,y,z] = azidip2plane(x0, y0, z0, strike, dip, RS, RD)
% Compute x/y/z coordinates of point on a plane with given strike and dip,
% at an along-strike distance RS, and a along-dip distance RD.
% Input units: angles in degrees, distances in meters
% 
% Transformation sketch in fear/proj/FIMoS/fig/trigonometry_sketch.JPG


dipdir = (strike-90)*pi/180;
dip    = dip        *pi/180;
strike = strike     *pi/180;

% Formulation for N/E/D
x = x0 + RS*cos(strike) + RD*cos(dip)*cos(dipdir);
y = y0 + RS*sin(strike) + RD*cos(dip)*sin(dipdir);
z = z0 + RD*sin(dip);

% Formulation for E/N/U
% x = x0 + RS*sin(strike) + RD*cos(dip)*sin(dipdir);
% y = y0 + RS*cos(strike) + RD*cos(dip)*cos(dipdir);
% z = z0                  - RD*sin(dip);


