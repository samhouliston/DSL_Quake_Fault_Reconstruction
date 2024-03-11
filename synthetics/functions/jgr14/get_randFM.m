function [strike,dip,rake] = get_randFM(nFM)

% Function for generating <nFM> random focal mechanism, after Kagan 2005
% All angles in degrees
%
% mam, 120623

strike = rand(nFM,1)*360;
dip    = acos(rand(nFM,1))*180/pi;
rake   = rand(nFM,1)*360-180;

strike = round(strike);
dip    = round(dip);
rake   = round(rake);

idx = find(strike==360);
if (~isempty(idx))
    strike(idx) = 0;
end