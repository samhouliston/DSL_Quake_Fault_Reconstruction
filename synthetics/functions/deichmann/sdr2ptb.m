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

x(abs(x)<1e-9) = 0; %mam
y(abs(y)<1e-9) = 0; %mam
z(abs(z)<1e-9) = 0; %mam

p = (sqrt(2)/2) .* (x - y);
t = (sqrt(2)/2) .* (x + y);

[paz,pdip] = dircos2azdip(p);
[taz,tdip] = dircos2azdip(t);
[baz,bdip] = dircos2azdip(z);

%Print formated:
%Str = sprintf('%03.0f/%02.0f %03.0f/%02.0f',paz,pdip,taz,tdip);
%disp(Str)