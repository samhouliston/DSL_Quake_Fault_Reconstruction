function [as,stress] = sample_rotatedEllipsexm1_stress_profile(nrnd,dtau0,dtaumax)


% NOT DONE



ang=0:0.01:2*pi; 

% Ellipse
a = 1;
b = 2;
psi = pi/6;
x = a*cos(ang)*cos(psi) - b*sin(ang)*sin(psi);
y = b*cos(psi)*sin(ang) + a*cos(ang)*sin(psi);
clf; hold on; grid on; box on; plot(x,y,'color','r');


azrnd = 2*pi*rand(nrnd,1); %max(azrnd); clf; histogram(azrnd)
xrnd  = a*cos(azrnd)*cos(psi) - b*sin(azrnd)*sin(psi);
yrnd  = b*sin(azrnd)*cos(psi) + a*cos(azrnd)*sin(psi);
rrnd  = sqrt(xrnd.^2+yrnd.^2);
for irnd = 1:nrnd
    
end

r = sqrt(x.^2+y.^2);
clf; hold on; grid on; box on; plot(ang,r,'color','r');




% rad=1;
% ang=0:0.01:2*pi; 
% xc=rad*cos(ang);
% yc=2*rad*sin(ang);
% 
% clf; hold on; axis equal
% set(gca,'xlim',[-2 2],'ylim',[-2 2])
% plot(xc,yc,'color','r');
% 
% rmin = 1;
% rmax = 2;
% r    = cos(2*ang)./(2*(rmax-rmin))+rmin+.5;
% clf; plot(ang,r)
% 
% 
% 
% xc=r.*cos(ang);
% yc=r.*sin(ang);
% plot(xc,yc,'color','r');