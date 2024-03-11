function fmpar = beachball_mod(strike1,dip1,rake1,lon0,lat0,radius,fignum,ptax,color)
%
% beachball(strike1,dip1,rake1,fignum,ptax,color)
%
% Plots double-couple faultplane solution with quadrants with downward 
% first motions shaded (lower-hemisphere, equal-area stereographic 
% projection with radius = 1.
% Input parameters are strike, dip and rake (degrees) of one of the nodal
% planes.
% Optional parameters are figure number (default 1), option to plot 
% P- and T-axes as red triangles (ptax = 'y' or 'n', default 'n') and
% fill-color (default 'k').
%
% beachball calls the following other functions:
%     sdr2sdr.m   calculates strike, dip and rake of conjugate plane
%     sdr2ptb.m   calculates azimuth and plunge of P-, T- and B-axis
%     pcircle.m   plots outer circle of stereo net
%     sd2xy.m     calculates x,y coordinates of nodal planes
%     arc2xy.m    calculates x,y coordinates of arc to close polygon
%     sd2norm     calculates vector components of strike, dip and normal
%     azdip2dircos.m  calculates vector components of B- and T-axis
%     equalarea.m calculates position of P- and T-axis on the beachball
%
% n. deichmann, sed-ethz, 2011/03/14.
%

if nargin < 9, color = 'k'; end
if nargin < 8, ptax = 'n'; end
if nargin < 7, fignum = 1; end

% Compute strike, dip and rake of conjugate nodal plane
[strike2,dip2,rake2] = sdr2sdr(strike1,dip1,rake1);

% Compute azimuth and plunge of P-, T- and B-axes
[paz,pdip,taz,tdip,baz,bdip] = sdr2ptb(strike1,dip1,rake1);

fmpar = sprintf('NP1:%3d/%2d/%4d  NP2:%3d/%2d/%4d  P-axis:%3d/%2d  T-axis:%3d/%2d',...
    strike1,dip1,rake1,round(strike2),round(dip2),round(rake2),...
    round(paz),round(pdip),round(taz),round(tdip));
fprintf('\n%s\n',fmpar)

figure(fignum);
hold on
pcircle_mod(lon0,lat0,radius,360,0,0,'k')
%
% Fault plane solution with three quadrants (pure normal or thrust)
%
if bdip < 0.5                              % horizontal B-axis
    if tdip < 45                           % pure normal fault
        [x,y] = sd2xy_mod(strike1,dip1,180,0,radius); % first nodal plane
        n1 = length(x);
        xx = x+lon0;
        yy = y+lat0;
        [x,y] = arc2xy_mod(strike1,strike1+180,radius);
        n2 = length(x);
        xx(n1+1:n1+n2) = x+lon0;
        yy(n1+1:n1+n2) = y+lat0;
        fill(xx,yy,color)
        [x,y] = sd2xy_mod(strike2,dip2,180,0,radius); % second nodal plane
        n1 = length(x);
        xx = x+lon0;
        yy = y+lat0;
        [x,y] = arc2xy_mod(strike2,strike2+180,radius);
        n2 = length(x);
        xx(n1+1:n1+n2) = x+lon0;
        yy(n1+1:n1+n2) = y+lat0;
        fill(xx,yy,color)
    else                                   % pure thrust fault
        [x,y] = sd2xy_mod(strike1,dip1,180,0,radius); % first nodal plane
        n1 = length(x);
        xx = x+lon0;
        yy = y+lat0;
        [x,y] = sd2xy_mod(strike2,dip2,180,0,radius); % second nodal plane
        n2 = length(x);
        xx(n1+1:n1+n2) = x+lon0;
        yy(n1+1:n1+n2) = y+lat0;
        fill(xx,yy,color)
    end 
    return
end
%
% Fault plane solution with four quadrants
%        
% Compute angles between strike of nodal planes and B-axis
b = azdip2dircos(baz,bdip);
[s1,d1,n1] = sd2norm(strike1,dip1);
deltab1    = acos(dot(s1,b)) * 180/pi;
[s2,d2,n2] = sd2norm(strike2,dip2);
deltab2    = acos(dot(s2,b)) * 180/pi;

% Compute angles between fault plane normals and T-axis
t = azdip2dircos(taz,tdip);
deltat1 = acos(dot(n1,t)) * 180/pi;
deltat2 = acos(dot(n2,t)) * 180/pi;

flag = 0;
if deltat1 > 90, flag = 1; end
if deltat2 > 90, flag = 1; end

if flag == 1     % angle between T-axis and a fault normal > 90 deg.
    %
    % first quadrant
    [x,y] = sd2xy_mod(strike1,dip1,0,deltab1,radius);
    n1 = length(x);
    xx = x+lon0;
    yy = y+lat0;
    [x,y] = sd2xy_mod(strike2,dip2,deltab2,0,radius);
    n2 = length(x);
    xx(n1+1:n1+n2) = x+lon0;
    yy(n1+1:n1+n2) = y+lat0;
    a2 = strike2;
    a1 = strike1;
    delta = a2 - a1;
    if delta >  180, a1 = a1 + 360; end
    if delta < -180, a2 = a2 + 360; end
    [x,y] = arc2xy_mod(a2,a1,radius);
    n2 = length(x);
    n1 = length(xx);
    xx(n1+1:n1+n2) = x+lon0;
    yy(n1+1:n1+n2) = y+lat0;
    fill(xx,yy,color)
    %
    % second quadrant
    [x,y] = sd2xy_mod(strike1,dip1,180,deltab1,radius);
    n1 = length(x);
    xx = x+lon0;
    yy = y+lat0;
    [x,y] = sd2xy_mod(strike2,dip2,deltab2,180,radius);
    n2 = length(x);
    xx(n1+1:n1+n2) = x+lon0;
    yy(n1+1:n1+n2) = y+lat0;
    a2 = strike2+180;
    if a2 > 360, a2 = a2 - 360; end
    a1 = strike1+180;
    if a1 > 360, a1 = a1 - 360; end
    delta = a2 - a1;
    if delta >  180, a1 = a1 + 360; end
    if delta < -180, a2 = a2 + 360; end
    [x,y] = arc2xy_mod(a2,a1,radius);
    n2 = length(x);
    n1 = length(xx);
    xx(n1+1:n1+n2) = x+lon0;
    yy(n1+1:n1+n2) = y+lat0;
    fill(xx,yy,color)
else               % angle between T-axis and both fault normals <= 90 deg.
    %
    % first quadrant
    [x,y] = sd2xy_mod(strike1,dip1,0,deltab1,radius);
    n1 = length(x);
    xx = x+lon0;
    yy = y+lat0;
    [x,y] = sd2xy_mod(strike2,dip2,deltab2,180,radius);
    n2 = length(x);
    xx(n1+1:n1+n2) = x+lon0;
    yy(n1+1:n1+n2) = y+lat0;
    a2 = strike2+180;
    if a2 > 360, a2 = a2 - 360; end
    a1 = strike1;
    delta = a2 - a1;
    if delta >  180, a1 = a1 + 360; end
    if delta < -180, a2 = a2 + 360; end
    [x,y] = arc2xy_mod(a2,a1,radius);
    n2 = length(x);
    n1 = length(xx);
    xx(n1+1:n1+n2) = x+lon0;
    yy(n1+1:n1+n2) = y+lat0;
    fill(xx,yy,color)
    %
    % second quadrant
    [x,y] = sd2xy_mod(strike1,dip1,180,deltab1,radius);
    n1 = length(x);
    xx = x+lon0;
    yy = y+lat0;
    [x,y] = sd2xy_mod(strike2,dip2,deltab2,0,radius);
    n2 = length(x);
    xx(n1+1:n1+n2) = x+lon0;
    yy(n1+1:n1+n2) = y+lat0;
    a2 = strike2;
    a1 = strike1+180;
    if a1 > 360, a1 = a1 - 360; end
    delta = a2 - a1;
    if delta >  180, a1 = a1 + 360; end
    if delta < -180, a2 = a2 + 360; end
    [x,y] = arc2xy_mod(a2,a1,radius);
    n2 = length(x);
    n1 = length(xx);
    xx(n1+1:n1+n2) = x+lon0;
    yy(n1+1:n1+n2) = y+lat0;
    fill(xx,yy,color)
end  
hold off
%
% plot position of P- and T-axis as red triangles
%
if ptax == 'y'
    hold on
    r = 0.04;
    [xx,yy]= equalarea(paz,90-pdip);
    pcircle(xx,yy,r,3,0,'n','r')       % P - empty triangle
    [xx,yy]= equalarea(taz,90-tdip);
    pcircle(xx,yy,r,3,0,'y','r')       % T - filled triangle
    hold off
end

function c = azdip2dircos(az,dip)
%
%     c = azdip2dircos(az,dip);
%
%     Converts azimuth and dip (degrees) of a vector into its direction 
%     cosines.
%     dip is actually the plunge of the vector.
%
%     The direction cosines are stored in array c.
%
%     Assume a right-handed coordinate system: 1 = N, 2 = E, 3 = Z(down).
%     The components of a unit-vector defined by azimuth, measured
%     counterclockwise from north, and dip, measured downward from the
%     horizontal, are then given by:
%     c(1) = cos(dip) * cos(azim)
%     c(2) = cos(dip) * sin(azim)
%     c(3) = sin(dip)
%
%     n. deichman, sed-ethz, 2011/03/13
%
rad = pi/180;
ra  = az*rad;
rd  = dip*rad;
c(1) = cos(rd) * cos(ra);
c(2) = cos(rd) * sin(ra);
c(3) = sin(rd);

function [az,dip] = dircos2azdip(c);
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
    theta = asin(abs(c(2)/codip));
else
    theta = acos(abs(c(1)/codip));
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

function [x,y]= equalarea(phi,theta)
%
% [x,y]= equalarea(phi,theta);
% 
% Given a point with position defined by its azimuth phi
% (counterclockwise from N) and its vertical angle theta (measured upward
% from the downgoing vertical axis) compute its x and y coordinates on 
% an equal area, lower hemisphere stereographic projection (Schmidt).
%
% phi and theta can also be vectors of multiple points.
%
% For points defined by their dip (measured downward from the horizontal)
% use [x,y]= equalarea(phi,90-dip);
%
% The stereonet is assumed to have a radius = 1.
%
% n. deichmann, sed-ethz, 2011/03/13.
%
alpha = phi;
beta  = theta;
i = find(beta > 90);
if length(i) > 0;
    beta(i) = 180 - beta(i);
    alpha(i) = alpha(i) -180;
end
radeg = pi/180;
z = sqrt(2) * sin(0.5*radeg.*(beta));
x = z .* sin(radeg.*alpha);
y = z .* cos(radeg.*alpha);

function [strike,dip] = norm2sd(n)
%
%     [strike,dip] = norm2sd(n);
%
%     Given the downward normal vector to a plane, calculate the 
%     strike and dip.
%     strike is measured counterclockwise from north and
%     dip is measured downward from the horizontal (degrees).
%
%     Assume a right-handed coordinate system: 1 = N, 2 = E, 3 = Z(down).
%
%     This function calls function dircos2azdip.
%
%     n. deichman, sed-ethz, 2011/03/13
%

%     s(1),s(2),s(3) are the dircos of strike;
%     d(1),d(2),d(3) are the dircos of dip;
%     n(1),n(2),n(3) are the dircos of the normal to the plane.
%     The vectors s, d, n form a right-handed coordinate system,
%     with n being the downward normal of the plane.
%
d(3) =  sqrt(1 - n(3)^2);  % Z-component of dip is perp to Z component of n
d(2) = -n(2) * n(3)/d(3);  % E-component of dip
d(1) = -n(1) * n(3)/d(3);  % N-component of dip
[azimuth,dip]   = dircos2azdip(d);  % dip of fault plane
s = cross(d,n);            % vector product: strike is perp to dip and norm
s(3) = 0;                  % force Z-component of strike to 0
[strike,plunge] = dircos2azdip(s);  % strike of fault plane

function [s,d,n] = sd2norm(strike,dip)
%
%     [s,d,n] = sd2norm(strike,dip);
%
%     Given strike and dip of a plane, calculate the direction cosines
%     of the vectors in direction of its strike, dip and normal.
%     Strike is measured counterclockwise from north and
%     dip is measured downward from the horizontal (degrees).
%
%     Assume a right-handed coordinate system: 1 = N, 2 = E, 3 = Z(down).
%
%     s(1),s(2),s(3) are the dircos of strike;
%     d(1),d(2),d(3) are the dircos of dip;
%     n(1),n(2),n(3) are the dircos of the normal to the plane.
%     The vectors s, d, n form a right-handed coordinate system,
%     with n being the downward normal of the plane.
%
%     This function calls function azdip2dircos.
%
%     n. deichman, sed-ethz, 2011/03/13
%
s = azdip2dircos(strike,0);      % vector of strike
d = azdip2dircos(strike+90,dip); % vector of dip
n = cross(s,d);                  % vector product: normal is perp to strike and dip

function [x,y] = sd2xy_mod(strike,dip,a1,a2,radius)
%
% [x,y] = sd2xy(strike,dip,a1,a2);
%
% Calculates x and y coordinates in an equal area stereographich projection
% (lower hemisphere, radius = 1) of a nodal plane segment with given strike
% and dip (degrees).
% The segment is given by the angles a1 and a2, measured in the plane
% downwards from the strike (degrees).
% For the complete nodal plane, set a1 = 0 and a2 = 180.
% 
% n. deichmann, sed-ethz, 2011/03/14.
%
radeg = pi/180;
theta1 = (90-a1)*radeg;
theta2 = (90-a2)*radeg;
sinazim = sin(strike*radeg);
cosazim = cos(strike*radeg);
%
% vertical nodal plane
%
if dip > 89.5                
    x = [0 0];
    y(1) = sign(theta1) * sqrt(1-cos(theta1));
    y(2) = sign(theta2) * sqrt(1-cos(theta2));
    xx = x;
    yy = y;
    x =  (xx.*cosazim + yy.*sinazim).*radius;
    y = (-xx.*sinazim + yy.*cosazim).*radius;
    return
end
%
% inclined nodal plane
%
sindip = sin(dip*radeg);
cosdip = cos(dip*radeg);
n = round(abs(theta1-theta2)/radeg);
if n < 2, n = 2; end
theta = linspace(theta1,theta2,n);
% calculate points for a nodal plane that strikes N-S
for i = 1:n
            sintheta = abs(sin(theta(i)));
            costheta = abs(cos(theta(i)));
            if sintheta < 0.0000001, sintheta = 0.0000001; end
            cotan = costheta/sintheta;
            z1 = 1 - costheta*sindip;
            z2 = cotan*cosdip;
            z2 = z2*z2;
            x(i) = sqrt(z1/(1 + 1/z2));
            y(i) = sign(theta(i)) * sqrt(z1/(1 + z2));
end
% rotate nodal plane to given strike
xx = x;
yy = y;
x =  (xx.*cosazim + yy.*sinazim).*radius;
y = (-xx.*sinazim + yy.*cosazim).*radius;

function [paz,pdip,taz,tdip,baz,bdip] = sdr2ptb(strike,dip,rake)
%
%     [paz,pdip,taz,tdip,baz,bdip] = sdr2ptb(strike,dip,rake);
%
%     Given strike, dip and rake of a fault plane calculate the slip vector.
%     strike is measured counterclockwise from north,
%     dip is measured downward from the horizontal and
%     rake is measured in the plane relative to the strike (degrees).
%     u is an array with the direction cosines of the slip vector.
%     u is also the normal to the conjugate fault plane (downward if
%     rake is negative and upward if rake is positive).
%
%     Assume a right-handed coordinate system: 1 = N, 2 = E, 3 = Z(down).
%
%     This function calls function dircos2azdip.
%
%     n. deichman, sed-ethz, 2011/03/13
%
radeg = pi/180;
sistrike = sin(strike*radeg);
costrike = cos(strike*radeg);
sidip    = sin(dip*radeg);
codip    = cos(dip*radeg);
sirake   = sin(rake*radeg);
corake   = cos(rake*radeg);

y(1) = -corake*costrike - sirake*codip*sistrike;
y(2) = -corake*sistrike + sirake*codip*costrike;
y(3) =  sidip*sirake;
x(1) =  sistrike*sidip;
x(2) = -costrike*sidip;
x(3) =  codip;
z(1) =  costrike*sirake - corake*codip*sistrike;
z(2) =  sistrike*sirake + corake*codip*costrike;
z(3) =  sidip*corake;
if z(3) < 0
    z = -z;
end

p = (sqrt(2)/2) .* (x - y);
t = (sqrt(2)/2) .* (x + y);

[paz,pdip] = dircos2azdip(p);
[taz,tdip] = dircos2azdip(t);
[baz,bdip] = dircos2azdip(z);

function [strike2,dip2,rake2] = sdr2sdr(strike1,dip1,rake1)
%
%     [strike2,dip2,rake2] = sdr2sdr(strike1,dip1,rake1);
%
%     Given strike, dip and rake of a fault plane (1) calculate 
%     strike, dip and rake of the conjugate fault plane (2).
%     strike is measured counterclockwise from north,
%     dip is measured downward from the horizontal and
%     rake is measured in the plane relative to the strike (degrees).
%
%     Assume a right-handed coordinate system: 1 = N, 2 = E, 3 = Z(down).
%
%     This function calls functions sd2norm, sdr2slip, norm2sd and
%     azdip2dircos.
%
%     n. deichman, sed-ethz, 2011/03/14
%

% case of horizontal conjugate plane
if dip1 > 89.5                    % dip1 vertical
    if abs(rake1) > 89.5          % rake1 vertical
        rake2 = rake1;       
        dip2 = 0;                 % dip2 horizontal
        strike2 = strike1 - 180;
        if strike2 < 0, strike2 = strike1 + 180; end
        return
    end
end

% all other cases
[s1,d1,n1] = sd2norm(strike1,dip1);  % vector components of strike, dip and normal
if rake1 < 0                         % normal faults
    n2 = sdr2slip(strike1,dip1,rake1);     % normal to plane2 = slip on plane1
else                                 % strike-slip and thrust faults
    n2 = sdr2slip(strike1,dip1,rake1-180); % normal must point downwards
end
[strike2,dip2] = norm2sd(n2);        % strike and dip of plane2
s2 = azdip2dircos(strike2,0);        % vector components of strike2
u2 = n1;                             % vector components of slip on plane2
rake2 = -acos(dot(u2,s2)) * 180/pi;  % scalar product gives angle between u2 and s2
if rake1 >= 0                        
    rake2 = rake2 + 180;             % if thrust mechanism, flip rake2
end

function u = sdr2slip(strike,dip,rake)
%
%     u = sdr2slip(strike,dip,rake);
%
%     Given strike, dip and rake of a fault plane calculate the slip vector.
%     strike is measured counterclockwise from north,
%     dip is measured downward from the horizontal and
%     rake is measured in the plane relative to the strike (degrees).
%     u is an array with the direction cosines of the slip vector.
%     u is also the normal to the conjugate fault plane (downward if
%     rake is negative and upward if rake is positive).
%
%     Assume a right-handed coordinate system: 1 = N, 2 = E, 3 = Z(down).
%
%     n. deichman, sed-ethz, 2011/03/13
%
radeg = pi/180;
sistrike = sin(strike*radeg);
costrike = cos(strike*radeg);
sidip    = sin(dip*radeg);
codip    = cos(dip*radeg);
sirake   = sin(rake*radeg);
corake   = cos(rake*radeg);
u(1) =  corake*costrike + sirake*codip*sistrike;
u(2) =  corake*sistrike - sirake*codip*costrike;
u(3) =  sidip*sirake;

function [x,y] = arc2xy_mod(azim1,azim2,radius)
%
% function [x,y] = arc2xy(azim1,azim2);
%
% Calculates x,y coordinates of an arc along a circle with radius = 1,
% centered at 0,0.
% The segment is delimited by the angles azim1 and azim2 measured
% clockwise from north (degrees) and the angle increment is 1 degree.
% It is used by function beachball.m to close the polygons of the shaded
% quadrants of the fault plane solution.
%
% n. deichmann, sed-ethz, 2011/03/14.
%
radeg = pi/180;
a1 = azim1 * radeg;
a2 = azim2 * radeg;
n = round(abs(a1-a2)/radeg);
if n < 2, n = 2; end
a = linspace(a1,a2,n);
x = radius*sin(a);
y = radius*cos(a);

function pcircle (xo,yo,radius,npoints,angle,ifill,color)
%
% Function PCIRCLE plots a circle as polygon.
%
% pcircle(xo,yo,radius,npoints,angle,ifill,color)
%
% xo,yo   = coordinates of center (def = 0,0)
% radius  = length of radius (def = 1)
% npoints = number of points on circle (def = circle)
% angle   = rotation of starting point in degrees (def = top)
% ifill   = fill circle ('n' --> no, 'y' --> yes) (def = no)
% color   = color as in function plot ('k','b','r',...) (def = black)
%
% author: n.deichmann, geophysics, eth-z, 2006/05/08.
% 
if nargin < 7, color = 'k'; end
if nargin < 6, ifill = 'n'; end
if nargin < 5, angle = 0;   end
if nargin < 4, npoints = 0; end
if nargin < 3, radius = 1;  end
if nargin < 2, yo = 0;      end
if nargin < 1, xo = 0;      end
%
resolution = 0.02;
if radius==0, radius = 1; end
if npoints == 0, npoints = round(radius * 2*pi/resolution); end
if npoints > 360, npoints = 360; end
n = npoints;
if n<3, error ('n < 3, no plot'), end
phi0 = pi/2 + angle*pi/180;
%
phi  = linspace(phi0,phi0+2*pi,n+1);
x = xo + radius * cos(phi);
y = yo + radius * sin(phi);
if ifill == 'y'
  fill (x,y,color)
else
  plot (x,y,color)
end  
axis equal off


