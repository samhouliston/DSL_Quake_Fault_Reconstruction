function [quakes,stress] = sample_Dieterich94_stress_profile(...
    neq, ...
    crack_radius, ...
    events_inside_crack, ...
    plotme)

% Compute stress profile from Dieterich 1994, eq. 20 and randomly sample
% from it, to simulate aftershock distributions as a function of distance 
% from the main shock hypocentre. Then distribute these distances on a 2D
% plane, by randomly sampling an azimuth.
% Men-Andrin Meier, 26/5/2023

                         
if nargin<4; plotme=false;
end
if nargin<3; events_inside_crack=false;
end

quakes.mname = 'Dieterich94';
quakes.crad  = crack_radius;




% 1. Compute stress profile
% Equation 20 of Dieterich 1994, JGR: "Simplified representation of stress 
% change from slip on a planar surface in a homogeneous elastic medium"
% Only defined for x>c
dx = 0.01;
x  = 0:dx:10*crack_radius;

dtauEq  = -2;
dtaumax = 5*abs(dtauEq);
dtau0   = -dtauEq*((1-crack_radius^3./x.^3).^(-1/2)-1);

% Equation only holds for x>c
% Set stress inside crack to 0:      no events will be sampled from within c
% Set stress inside crack to taumax: numerous events sampled from within c
dtau                 = dtau0;
if events_inside_crack; dtau(x<crack_radius) = dtaumax;
else                    dtau(x<crack_radius) = 0;
end
dtau(dtau>dtaumax) = dtaumax;




% 2. Use stress profile as pdf to sample random main shock / aftershock 
% distances  

% Compute cdf, and use uniformly random distribution [0, 1] to sample on
% y-axis
dtaui  = cumsum(dtau)*dx;
dtauin = dtaui./(sum(dtau)*dx); %clf; plot(x,dtauin)

frnd = rand(neq,1);
ri = zeros(neq,1);
for irnd = 1:neq
    [~,idx]  = min(abs(dtauin-frnd(irnd)));
    ri(irnd) = x(idx);
end

% Random position on 2D plane with uniformly random azimuth
azi = 2*pi*rand(neq,1); %max(azrnd); clf; histogram(azrnd)
xi  = ri.*cos(azi);
yi  = ri.*sin(azi);
%zi  = d3*randn(neq,1);
%clf; plot(xrnd,yrnd,'.k')

stress.x    = x;
stress.dtau = dtau;
quakes.ri   = ri;
quakes.xi   = xi;
quakes.yi   = yi;
%quakes.zi   = zi;


if plotme

    figure(1); clf; 
    subplot(1,2,1); hold on; grid on; box on;
    %plot(x,dtau0,'-k','lineWidth',1)
    plot(x,dtau ,'-r','lineWidth',2)
    set(gca,'yscale','log')
    histogram(ri, 0:.05:max(x))
    
    subplot(1,2,2); hold on; grid on; box on; axis equal
    plot(quakes.xi, quakes.yi, '.k')
    xlabel('x')
    ylabel('y')
end