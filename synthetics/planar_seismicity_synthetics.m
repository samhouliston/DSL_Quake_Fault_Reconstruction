clear all
addpath(genpath(pwd))

% This is a script for generating synthetic seismicity clouds around planar
% features with a given strike, dip and size. The script uses a Cartesian 
% North/East/Down coordinate system with metre units.

set(groot,'DefaultFigureWindowStyle','docked')
set(groot,'defaultAxesFontName','Arial')
set(groot,'defaultAxesFontSize',10)



%% Specify plane orientation
stk0 = 15; 
dip0 =  80;

% Specify hypocentral coordiantes, magnitude & stressdrop
n0  = 0; 
e0  = 0;
d0  = 0;
mag = 1; % Source size
stressdrop = 1e6;  % [Pa]

plane_colour = [1 0 1];
view_angle   = [-115 35];





%% Compute rectangular finite source patch
finsrc = get_rectangular_slippatch_from_FM(...
    n0, e0, d0, ...
    stk0, ...
    dip0, ...
    mag, ...
    stressdrop);



% Plot
figure(2001); clf; hold on; grid on; box on; 
plot_rectangular_slippatch(finsrc, plane_colour)

xlabel('Northing [m]')
ylabel('Easting [m]')
zlabel('Down [m]')

lim = get_cubic_domain_XYZ(...
    n0, ...
    e0, ...
    d0, ...
    1.5*finsrc.length);

set(gca, ...
    'yDir','reverse', ...
    'zDir','reverse', ... % used for NEDown
    'view', view_angle, ...
    'xlim', lim.x, ...
    'ylim', lim.y, ...
    'zlim', lim.z)




%% Create synthetic seismicity catalogue

% Use Dieterich 1994 stress decay model to sample spatial event densities
% Then start with a reference FM, create <neq> replicas and perturb them with
% a random rotation angle around a fixed or random rotation axis
neq     = 1e3;           % Number of quakes
drot    = 15;            % Rotation angle variability
cradius = finsrc.length; % Crack radius, definse cluster size
dx      = cradius/25;    % Stdev of scatter around perfect plane

% Generate hypocentres
cat = get_seismicity_sample_from_Dieterich94( ...
            neq, ...
            n0, e0, d0, ...
            stk0, dip0, ...
            cradius, ...
            dx); 

% Plot 
plot3(cat.n, cat.e, cat.d, '.', 'color', [.2 .2 .2])

% Select plotting box
plt.lim = set_bounding_box([cat.n, cat.e, cat.d], .4, 99.6);






%% Infer cluster orientation with Principal Component Analysis
evd = get_eigenval_decomp_NED(cat);
plot_eigenvectors(evd, n0, e0, d0)


% Infer plane orientation from eigenvectors
% 3rd eigenvector is normal vector to inferred plane
[stkHat,dipHat] = norm2sd(evd.V(:,3));

[stk0, stkHat, stk0-stkHat; ...
 dip0, dipHat, dip0-dipHat]

tstring1 = sprintf('true     fault orientation: strike=%i°, dip=%i°', stk0, dip0);
tstring2 = sprintf('inferred fault orientation: strike=%i°, dip=%i°', round(stkHat), round(dipHat));
title({tstring1; tstring2});



%% Test strike & dip estimates with gridsearch
if false; compare_true_and_inferred_strike_dip;
end