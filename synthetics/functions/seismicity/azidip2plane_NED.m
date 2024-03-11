function [n,e,d] = azidip2plane_NED(n0, e0, d0, strike, dip, RS, RD)
% Compute n/e/d coordinates of point on a plane with given strike and dip,
% at an along-strike distance RS, and a along-dip distance RD.
% Input units: angles in degrees, distances in meters
% 
% Transformation sketch in fear/proj/FIMoS/fig/trigonometry_sketch.JPG


dipdir = (strike-90)*pi/180;
dip    = dip        *pi/180;
strike = strike     *pi/180;

% Formulation for N/E/D
n = n0 + RS*cos(strike) + RD*cos(dip)*cos(dipdir);
e = e0 + RS*sin(strike) + RD*cos(dip)*sin(dipdir);
d = d0                  - RD*sin(dip);

% Formulation for E/N/U
% n = n0 + RS*sin(strike) + RD*cos(dip)*sin(dipdir);
% e = e0 + RS*cos(strike) + RD*cos(dip)*cos(dipdir);
% d = d0                  - RD*sin(dip);


