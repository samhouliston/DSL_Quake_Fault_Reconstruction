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