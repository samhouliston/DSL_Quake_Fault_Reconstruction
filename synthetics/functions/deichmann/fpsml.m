function fpsml(finam,strike1,dip1,rake1,txt)
%
% fpsml(finam,strike1,dip1,rake1,txt)
%
% Function to plot first motion polarities and focal mechanisms on a
% stereo-net (lower hemisphere, equal area, with radius = 1).
% Input:
% finam:  filename of nonlinloc output (.MANUPDE format)
%         (if not available, give ' ' (one blank character) instead)
% strike1,dip1,rake1: strike, dip and rake of one nodal plane
%         (if not given, only first motions will be plotted)
% txt:    'n' = no text on plot; 'y' = station names and parameters on plot
%         (default is 'y').
%
% fpsml computes strike,dip and rake of the second nodal plane as well as 
% P- and T-axes for the focal mechanism defined by strike1, dip1 and rake1. 
% If txt = 'y', the positions of the P- and T-axis are plotted as red 
% triangles.
%
% If both an input file with first motions as well as strike, dip and rake
% are given, the output will be written to file fpsml.txt and a summary
% line will be appended to file fpsml.log.
%
% fpsml calls the following other functions:
%    pdeloc.m and pdestn.m to read the input data from the .MANUPDE file
%    pcircle.m   to draw circle of stereo-net and first motion polarities
%    equalarea.m to calculate x and y coordinates of points on stereo-net
%    sdr2sdr.m   to compute strike, dip and rake of second nodal plane
%    sdr2ptb.m   to compute azimuth and plunge of P-, T- and B-axes
%    sd2xy.m     to draw the nodal planes on the stereo net
%    faultregime.m to determine fault regime according to World Stress Map
%    thetaphi.m  to compute take-off angles relative to normal and strike
%    radpat.m    to compute radiation coefficients for P, SV and SH waves
%
% n. deichmann, sed-ethz, 2011/03/15.
%
clf
pcircle                        % plot outer circle of stereo-net
%
if nargin < 5, txt = 'y'; end   % default: plot text and station names
%
nstn = 0;
if finam(1) ~= ' '
    %
    % Read event location from input file
    [yr,mo,dy,hr,mi,sec,lat,lon,depth,mag,dm,rms,gap,no]   = pdeloc(finam);
    %
    % Read stations, polarities and take-off angles from input file
    [stn,pha,pol,atm,ats,res,amp,smag,umag,dist,azim,angl] = pdestn(finam);
    %
    % Compute x and y coordinates on stereo-net for all stations
    [x,y]= equalarea(azim,angl); 
    %
    % Plot available polarities on stereo-net
    r = 0.04;
    if txt == 'y', r=0.03; end
    nstn = length(x);
    hold on
    for i = 1:nstn
        if pol(i,2) ~= ' '
          if pol(i,2) == 'U', pcircle(x(i),y(i),r,36,0,'y','k'), end  % filled circles
          if pol(i,2) == 'D', pcircle(x(i),y(i),r,36,0,'n','k'), end  % empty circles
          if (txt == 'y')
            angle = atan2(y(i),x(i))*180/pi;
            text (x(i)+2*x(i)*r,y(i)+2*y(i)*r,stn(i,1:5),'Rotation',[angle])
          end
        end
    end
    %
    % Add text with event parameters
    evpar = sprintf('%4d/%2d/%2d %2d:%2d:%4.1f  %6.3f  %6.3f  %4.1f  %3.1f',...
        yr,mo,dy,hr,mi,sec,lat,lon,depth,mag);
    fprintf('\n%s\n',evpar);
    if (txt == 'y')
        text(-1.0,1.2,evpar)
        plot (0,0,'k.')
    end
    hold off
    axis tight
    %
    % Position figure in lower left corner of screen
    set(gcf,'Position',[10 40 560 580])
    set(gcf,'PaperPositionMode','auto')
end  
%
if nargin < 4, return, end   % strike, dip and rake are not given
%
% -----------------------------------------------------------------------
%
% Compute strike, dip and rake of conjugate nodal plane
[strike2,dip2,rake2] = sdr2sdr(strike1,dip1,rake1);
%
% Plot the two nodal planes
hold on
[x,y] = sd2xy(strike1,dip1,0,180);
plot(x,y,'k')
[x,y] = sd2xy(strike2,dip2,0,180);
plot(x,y,'k')
hold off
%
% Compute azimuth and plunge of P- and T-axes
[paz,pdip,taz,tdip,baz,bdip] = sdr2ptb(strike1,dip1,rake1);
%
% plot position of P- and T-axes as red triangles if txt='y'
if (txt == 'y')
    hold on
    r = 0.04;
    [xx,yy]= equalarea(paz,90-pdip);
    pcircle(xx,yy,r,3,0,'n','r')       % P - empty triangle
    [xx,yy]= equalarea(taz,90-tdip);
    pcircle(xx,yy,r,3,0,'y','r')       % T - filled triangle
    hold off
end
%
% Add text with focal mechanism parameters
fmpar = sprintf('NP1:%3d/%2d/%4d  NP2:%3d/%2d/%4d  P-axis:%3d/%2d  T-axis:%3d/%2d',...
    strike1,dip1,rake1,round(strike2),round(dip2),round(rake2),...
    round(paz),round(pdip),round(taz),round(tdip));
fprintf('\n%s\n',fmpar);
if (txt == 'y')
    hold on
    text(-1.0,-1.2,fmpar)
    hold off
end
axis tight
%
% Position figure in lower left corner of screen
set(gcf,'Position',[10 40 560 580])
set(gcf,'PaperPositionMode','auto')
%
if nstn == 0, return, end
%
% ------------------------------------------------------------------------
%
% If input file and strike, dip and rake are given, write results to
% output file fpsml.txt
fid1 = fopen('fpsml.txt','w+t');
fprintf(fid1,'\n%s\n',evpar);
fprintf(fid1,'\n%s\n',fmpar);
%
% Determine faulting regime according to world stress map criteria
regime = faultregime(pdip,tdip,bdip);
fprintf('\nFaulting regime (World Stress Map, JGR, 97, 11711, 1992): %s\n',regime);
fprintf(fid1,'\nFaulting regime (World Stress Map, JGR, 97, 11711, 1992): %s\n',regime);
%
% Compute angles theta and phi and radiation coefficients
% fprintf('\nSTAT     DIST  AZIM  ANGL  THETA   PHI      RP     RSV     RSH\n');
fprintf(fid1,'\nSTAT     DIST  AZIM  ANGL  THETA   PHI      RP     RSV     RSH\n');
for i = 1:nstn
    if pol(i,2)~=' '
        [theta,phi] = thetaphi(strike1,dip1,azim(i),angl(i));
        [rp,rsv,rsh] = radpat(strike1,dip1,rake1,azim(i),angl(i));
%        fprintf('%s %s %4d  %4d  %4d   %4d  %4d   %6.3f  %6.3f  %6.3f\n',...
%            stn(i,1:6),pol(i,2),round(dist(i)),round(azim(i)),round(angl(i)),...
%            round(theta),round(phi),rp,rsv,rsh);
        fprintf(fid1,'%s %s %4d  %4d  %4d   %4d  %4d   %6.3f  %6.3f  %6.3f\n',...
            stn(i,1:6),pol(i,2),round(dist(i)),round(azim(i)),round(angl(i)),...
            round(theta),round(phi),rp,rsv,rsh);
    end
end
fclose(fid1);
%
% Append event and focal mechanism parameters to fpsml.log
fid2 = fopen('fpsml.log','a+t');
fprintf(fid2,'%4d/%2d/%2d %2d:%2d %6.3f %6.3f %4.1f %3.1f %3d/%2d/%4d %3d/%2d/%4d %3d/%2d %3d/%2d\n',...
    yr,mo,dy,hr,mi,lat,lon,depth,mag,...
    strike1,dip1,rake1,round(strike2),round(dip2),round(rake2),...
    round(paz),round(pdip),round(taz),round(tdip));
fclose(fid2);
%
% ------------------------------------------------------------------------

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

function regime = faultregime(pdip,tdip,bdip)
%
% regime = faultregime(pdip,tdip,bdip);
%
% Given plunge of P-, T- and B-axes, determine faulting regime according 
% to World Stress Map criteria (JGR, 97, 11711, 1992).
%
% regime is either U, NF, TF, NS, TS or SS.
%
% n. deichmann, sed-ethz, 2011/03/14.
%
regime = 'U ';
if round(pdip) >= 52
    if round(tdip) <= 35, regime = 'NF'; end
end
if round(tdip) >= 52
    if round(pdip) <= 35, regime = 'TF'; end
end
if round(pdip) >= 40
    if round(pdip) < 52
        if round(tdip) <= 20, regime = 'NS'; end
    end
end
if round(tdip) >= 40
    if round(tdip) < 52
        if round(pdip) <= 20, regime = 'TS'; end
    end
end
if round(pdip) < 40
    if round(bdip) >= 45
        if round(tdip) <= 20, regime = 'SS'; end
    end
end
if round(tdip) < 40
    if round(bdip) >= 45
        if round(pdip) <= 20, regime = 'SS'; end
    end
end

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

function [x,y] = sd2xy(strike,dip,a1,a2)
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
    x =  xx.*cosazim + yy.*sinazim;
    y = -xx.*sinazim + yy.*cosazim;
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
x =  xx.*cosazim + yy.*sinazim;
y = -xx.*sinazim + yy.*cosazim;

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

function [rp,rsv,rsh] = radpat (strike,dip,rake,azimuth,takeoff)
%
% [rp,rsv,rsh] = radpat(strike,dip,rake,azimuth,takeoff);
%
% Function radpat calculates radiation coefficients of a 
% double couple source for P, SV and SH waves 
% (Aki & Richards, 2002, eqns. 4.89-4.91).
% 
% (Modified from Fortran program fpsml.f, 2007/04/30, n.d., sed-ethz)
%

      rad   = pi/180;
      delta = rad * dip;
      lamda = rad * rake;
      ixi   = rad * takeoff;
      dphi  = rad * (azimuth - strike);

      sinlamda  = sin(lamda);
      coslamda  = cos(lamda);
      sindelta  = sin(delta);
      cosdelta  = cos(delta);
      sindphi   = sin(dphi);
      cosdphi   = cos(dphi);
      sinixi    = sin(ixi);
      cosixi    = cos(ixi);
      sin2delta = sin(delta+delta);
      cos2delta = cos(delta+delta);
      sin2dphi  = sin(dphi+dphi);
      cos2dphi  = cos(dphi+dphi);
      sin2ixi   = sin(ixi+ixi);
      cos2ixi   = cos(ixi+ixi);

rp  = coslamda*sindelta*sinixi*sinixi*sin2dphi - ...
    coslamda*cosdelta*sin2ixi*cosdphi + ...
    sinlamda*sin2delta*(cosixi*cosixi - sinixi*sinixi*sindphi*sindphi) + ...
    sinlamda*cos2delta*sin2ixi*sindphi ;

rsv = sinlamda*cos2delta*cos2ixi*sindphi - ...
    coslamda*cosdelta*cos2ixi*cosdphi + ...
    0.5*coslamda*sindelta*sin2ixi*sin2dphi - ...
    0.5*sinlamda*sin2delta*sin2ixi*(1 + sindphi*sindphi) ;

rsh = coslamda*cosdelta*cosixi*sindphi + ...
    coslamda*sindelta*sinixi*cos2dphi + ...
    sinlamda*cos2delta*cosixi*cosdphi - ...
    0.5*sinlamda*sin2delta*sinixi*sin2dphi ;

function [theta,phi] = thetaphi(strike,dip,azim,angl)
%
% [theta,phi] = thetaphi(strike,dip,azim,angl);
%
% Given strike and dip of the active fault plane and the take-off angles
% azim and angl of a ray leaving the source, calculate the angle theta
% between the ray and the fault normal and the angle phi between the
% projection of the ray onto the fault and the strike of the fault.
%
% For explanations of the angles theta and phi see,
% Courboulex, Deichmann and Gariel, GJI, 139, 152-160.
%
% thetaphi calls functions sd2norm.m and azdip2dircos.m
%
% Modified from fortran subroutine thetaphi.f in program fpsml.f
% n. deichmann, sed-ethz, 2011/03/14
%
[s,d,n] = sd2norm(strike,dip);   % vector components of strike, dip and normal of fault plane
a = azdip2dircos(azim,90-angl);  % vector components of ray at the source in N,E,Z system
%
% rotate ray components into fault-based coordinate system (scalar products)
%
b(1) = a(1) * s(1) + a(2) * s(2) + a(3) * s(3);
b(2) = a(1) * d(1) + a(2) * d(2) + a(3) * d(3);
b(3) = a(1) * n(1) + a(2) * n(2) + a(3) * n(3);
%
% calculate angle between ray and normal to fault plane
%
theta = acos(abs(b(3))) * 180/pi;
%
% calculate phi
%
if theta < 0.5
    phi = 0;
    return
end
phi = atan2(b(2),b(1)) * 180/pi;

function [yr,mo,dy,hr,mi,sec,lat,lon,depth,mag,dm,rms,gap,no] = pdeloc(finam)

fid = fopen(finam);
ctxt = fgetl(fid);
if length(ctxt)<7, ctxt(1:7) = '       ';, end
while strcmp(ctxt(4:7),'DATE') == 0
    ctxt = fgetl(fid);
    if length(ctxt)<7, ctxt(1:7) = '       ';, end
end
ctxt = fgetl(fid);
disp (ctxt)
yr    = str2double(ctxt(2:5));
mo    = str2double(ctxt(6:7)); 
dy    = str2double(ctxt(8:9));
hr    = str2double(ctxt(11:12));
mi    = str2double(ctxt(14:15));
sec   = str2double(ctxt(17:20));
lat   = str2double(ctxt(22:27));
lon   = str2double(ctxt(30:36));
depth = str2double(ctxt(39:42));
dm    = str2double(ctxt(48:53));
mag   = str2double(ctxt(55:57));
rms   = str2double(ctxt(59:63));
gap   = str2double(ctxt(65:67));
no    = str2double(ctxt(69:71));
fclose(fid);

function [stn,pha,pol,atm,ats,res,amp,smag,umag,dist,azim,angl] = pdestn(finam)

fid = fopen(finam);
ctxt = fgetl(fid);
if length(ctxt)<3, ctxt(1:3) = '   ';, end
while strcmp(ctxt(1:3),'STA') == 0
    ctxt = fgetl(fid);
    if length(ctxt)<3, ctxt(1:3) = '   ';, end
end
i = 0;
while feof(fid)==0
  ctxt = fgetl(fid);
%  disp (ctxt)
  i = i + 1;
  stn(i,1:6) = ctxt(1:6);
  pha(i,1:7) = ctxt(8:14);
  pol(i,1:2) = ctxt(16:17);
  atm(i) = str2double(ctxt(21:22));
  ats(i) = str2double(ctxt(23:28));
  res(i) = str2double(ctxt(36:40));
  amp(i) = str2double(ctxt(51:58));
  if ctxt(69:69)==' '
      smag(i) = -9.9;
      umag(i,1) = '?';
  else
      smag(i) = str2double(ctxt(66:69));
      umag(i,1) = ctxt(70);
  end
  dist(i) = str2double(ctxt(72:78));
  azim(i) = str2double(ctxt(80:84));
  angl(i) = str2double(ctxt(86:90));
end
fclose(fid);