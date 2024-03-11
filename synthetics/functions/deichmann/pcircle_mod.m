function pcircle_mod(xo,yo,radius,npoints,angle,ifill,color)
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
%axis equal off