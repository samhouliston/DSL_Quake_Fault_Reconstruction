function [az,dip] = dircos2azdip(c)
%
%     [az,dip] = dircos2azdip(c);
%
%     Converts direction cosines of a vector into azimuth and dip (degrees).
%     dip is actually the plunge of the vector.
%
%     The direction cosines are stored in array c.
%     Assume a right-handed coordinate system: 1 = N, 2 = E, 3 = Z(down).
%     The components of a unit-vector defined by azimuth, measured
%     counterclockwise from north, and dip, measured downward from the
%     horizontal, are then given by:
%     c(1) = cos(dip) * cos(azim)
%     c(2) = cos(dip) * sin(azim)
%     c(3) = sin(dip)
%
%     if c(3) > 0.99999 then dip = 90 and az = 0!
%
%     n. deichman, sed-ethz, 2011/03/13
%
if c(3) > 0.99999
    dip = 90;
    az  =  0;
    return
end

delta = asin(c(3));
codip = cos(delta);
if abs(c(2)) > abs(c(1))
    y = abs(c(2)/codip); %mam - y was sometimes out of range [-1,1]
    if y>1; y=1; end     %mam   and returned complex number
    theta = asin(y);
else
    y = abs(c(1)/codip); %mam
    if y>1; y=1; end     %mam
    theta = acos(y);
end

dip = delta * 180/pi;
az  = theta * 180/pi;

if c(1) < 0
    if c(2) < 0
        az = 180 + az;
    else
        az = 180 - az;
    end
elseif c(2) < 0
    az = 360-az;
end
if dip < 0
    dip = -dip;
    az = az + 180;
    if az > 360
        az = az - 360;
    end
end
