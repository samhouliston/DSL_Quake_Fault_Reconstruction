function [as,stress] = sample_centerxm1_stress_profile(nrnd,dtau0,dtaumax)

plotme   = 1;
as.mname = 'rm1_center_wTaumax';
% Equation 20 of Dieterich 1994, JGR: "Simplified representation of stress 
% change from slip on a planar surface in a homogeneous elastic medium"
%c  = 1;         % Crack radius
dx   = 0.01;
x    = 0:dx:5;
dtau = dtau0./x.^3;
dtau(dtau>dtaumax) = dtaumax;  %clf; plot(x,dtau)

% Use dtau as pdf for aftershock radius; sample random locations 
dtaui  = cumsum(dtau)*dx;
dtauin = dtaui./(sum(dtau)*dx); %clf; plot(x,dtauin)

frnd = rand(nrnd,1);
rrnd = zeros(nrnd,1);
for irnd = 1:nrnd
    [~,idx] = min(abs(dtauin-frnd(irnd)));
    rrnd(irnd) = x(idx);
end

% Random position on 2D plane with uniformly random azimuth
azrnd = 2*pi*rand(nrnd,1); %max(azrnd); clf; histogram(azrnd)
xrnd  = rrnd.*cos(azrnd);
yrnd  = rrnd.*sin(azrnd);
zrnd  = .2*randn(nrnd,1); %clf; plot(xrnd,yrnd,'.k')

if plotme
    figure(333); clf; hold on; grid on; box on;
    plot(x,dtau0,'-r','lineWidth',1)
    plot(x,dtau ,'-k','lineWidth',2)
    set(gca,'yscale','log')
    
    histogram(xrnd,0:.05:max(x))
end

stress.x    = x;
stress.dtau = dtau;
as.dx       = xrnd;
as.dy       = yrnd;
as.dz       = zrnd;