function [n1,n2,n3,s1,s2,s3,b1,b2,b3] = fm2nsb_vectors(strike,dip,rake)
% Computes normal, slip and b-vectors for one of the two planes of a FM
% (strike,dip,rake in degrees)

deg2rad = pi/180;

strike = strike*deg2rad;
dip    = dip*deg2rad;
rake   = rake*deg2rad;

% % Stein and Wysession CS
% % normal vector
% n1 = -sin(dip)*sin(strike);
% n2 = -sin(dip)*cos(strike);
% n3 = -cos(dip);
% 
% % slip vector
% s1 = cos(rake)*cos(strike)+sin(rake)*cos(dip)*sin(strike);
% s2 = -cos(rake)*sin(strike)+sin(rake)*cos(dip)*cos(strike);
% s3 = sin(rake)*sin(dip);

% % focal-codes CS
% normal vector
n1 = -sin(dip)*sin(strike);
n2 =  sin(dip)*cos(strike);
n3 = -cos(dip);

% slip vector
s1 =  cos(rake)*cos(strike)+sin(rake)*cos(dip)*sin(strike);
s2 =  cos(rake)*sin(strike)-sin(rake)*cos(dip)*cos(strike);
s3 = -sin(rake)*sin(dip);

% % Copied from focal-codes
% anx=-sin(wdip)*sin(wstrik);
% any=sin(wdip)*cos(wstrik); % here I have a minus
% anz=-cos(wdip);
% 
% dx=cos(wrake)*cos(wstrik)+cos(wdip)*sin(wrake)*sin(wstrik); 
% dy=cos(wrake)*sin(wstrik)-cos(wdip)*sin(wrake)*cos(wstrik); % here I have both signs different
% dz=-sin(wdip)*sin(wrake); % here I have no minus


% b-vector, n x s = b
b1 = n2*s3 - n3*s2;
b2 = n3*s1 - n1*s3;
b3 = n1*s2 - n2*s1;

% If one of the components is close to 0, set it to zero
if ((n1 < 1E-5) && (n1 > -1E-5))
    n1 = 0;end
if ((n2 < 1E-5) && (n2 > -1E-5))
    n2 = 0;end
if ((n3 < 1e-05) && (n3 > -1E-5))
    n3 = 0;end
if ((s1 < 1E-5) && (s1 > -1E-5))
    s1 = 0;end
if ((s2 < 1E-5) && (s2 > -1E-5))
    s2 = 0;end
if ((s3 < 1e-05) && (s3 > -1E-5))
    s3 = 0;end
if ((b1 < 1E-5) && (b1 > -1E-5))
    b1 = 0;end
if ((b2 < 1E-5) && (b2 > -1E-5))
    b2 = 0;end
if ((b3 < 1e-05) && (b3 > -1E-5))
    b3 = 0;end