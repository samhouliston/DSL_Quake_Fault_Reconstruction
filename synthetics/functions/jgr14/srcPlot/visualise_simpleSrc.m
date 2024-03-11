function [Slon,Slat] = visualise_simpleSrc(ref_lat,ref_lon,ref_z,hyC_lat,hyC_lon,hyC_z,strike,s_length)
% computes coordinates of the 3 other corners of a source fault from hypocenter
% and reference point 3D coordinates

% UNITS: coordinates [deg], depths and s_length [km], strike [deg]
addpath(genpath('../../matlab/functions/src/'))

% Convert km to degrees ---------------------------------------------------
km2lat              = 1/111.2;
[one_deg_lon_in_km] = lon2km_smallCircle(1,hyC_lat);
km2lon              = 1/one_deg_lon_in_km;
% -------------------------------------------------------------------------

Neq      = 1;
Slat     = zeros(4,Neq);
Slon     = zeros(4,Neq);
% --> For both latitude and longitude, save one line per source as:
% [1,2]: upper start- and end-corner
% [3,4]: lower start- and end-corner
% [5,6]: start- and end-corner for map intersection

for i=1:1:Neq
    
    % Compute upper boundary coordinates of fault patch
    strike = pi()/180*strike;
    slat   = ref_lat;
    slon   = ref_lon;
    elat   = slat + slength*cos(sstrike)*km2lat;
    elon   = slon + slength*sin(sstrike)*km2lon;
    
    % Upper start- and end-corner
    Slat(1,i) = slat;             % Write pairs of point into matrices of
    Slat(2,i) = elat;             % dimensions 2 x Neq to combine coordinates
    Slon(1,i) = slon;             % with lines pairwise
    Slon(2,i) = elon;
    
    % Compute lower boundary coordinates of fault patch
    ds1 = hyC_lat - slat;
    ds2 = hyC_lon - slon;
    ds3 = hyC_lat - elat;
    ds4 = hyC_lon - elon;
    
    slat_low = elat + 2*ds3;
    slon_low = elon + 2*ds4;
    elat_low = slat + 2*ds1;
    elon_low = slon + 2*ds2;
    
    % Lower start- and end-corner
    Slat(3,i) = slat_low;             % Write pairs of point into matrices of
    Slat(4,i) = elat_low;             % dimensions 2 x Neq to combine coordinates
    Slon(3,i) = slon_low;             % with lines pairwise
    Slon(4,i) = elon_low;
    
end
