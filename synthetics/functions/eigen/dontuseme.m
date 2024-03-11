clear all
addpath(genpath(pwd))

% This script is for playing with beachballs. That is, to experiment with
% various operations and visualisation of earthquake focal mechanisms data.
% The script uses a North/East/Down coordinate system
% 
% Men-Andrin Meier, last updates 12/05/2023

% Includes functions from 
% . Nicholas Deichmanns toolbox (in fpsml.m)
% . Linus Villiger; from Villiger et al., 2021, GJI
% . Gregorz Kwiatek & Patricia Martinez-Garzon

% Open task list
% . triangle plots
% . seismicity clusters
% . PCA
% ...

set(groot,'DefaultFigureWindowStyle','docked')
set(groot,'defaultAxesFontName','Arial')
set(groot,'defaultAxesFontSize',10)



%% Specify focal mechanism
stk = 40; 
dip = 55;
rak = 60;

% Specify hypocentral coordiantes, magnitude & stressdrop
n0  = 0; 
e0  = 0;
d0  = 0;
mag = 1; % Source size
stressdrop = 1e6;  % [Pa]

% Plot properties
view_angle = [-100 75];



%% Compute all relevant vectors: slip, normals, P- and T-axes
fm = get_all_FM_vectors(stk,dip,rak);

% Check
if false; check_normal_vector_definition
end



%% Plot 3D beachball and finite source 
[hf1, fm, finsrc] = plot_3D_beachball_and_finite_source(...
    stk, dip, rak, ...
    n0, e0, d0, ...
    mag, ...
    stressdrop, ...
    view_angle);

if false; make_rotating_mp4;
          %make_rotating_mp4('~/fig/proj/beachballs/new/example_beachhball.mp4')
end


%% Lower hemisphere 2D visualisations
hf2 = plot_2D_FM_visualisations(fm);



 %% Create synthetic set of FMs
% Starting with a reference FM, create <nfm> replicas and perturb them with
% a random rotation angle around a fixed or random rotation axis
%stk0 = 30;  % Reference FM
%dip0 = 55;
%rak0 = 90;

nfm  = 400;  % Number of FMs
drot = 15;   % Rotation angle variability

FMs = get_FMs_from_rotation_around_random_axis(nfm, stk, dip, rak, drot);
%FMs = get_FMs_from_rotation_around_fixed_axis (nfm, stk, dip, rak, drot, [1 0 1]);

hf3 = plot_stereonet_for_set_of_FMs(FMs);

% Add reference mechanism
fm           = get_all_FM_vectors(stk, dip, rak);
[azi1, dip1] = norm_ned2azdip(fm.n1); 
[azi2, dip2] = norm_ned2azdip(fm.n2); 
[X1, Y1]     = stereoplot_point_cords(azi1, dip1);
[X2, Y2]     = stereoplot_point_cords(azi2, dip2);
h0 = plot([X1 X2], [Y1 Y2], 's', ...
    'markerFaceColor', 'c', ...
    'markerEdgeColor', [.4 .4 .4], ...
    'DisplayName', 'n_{ref. FM}');








%% Create synthetic seismicity cluster
neq     = 100;           % Number of earthquakes
cradius = finsrc.length; % Crack radius
dx      = cradius/20;    % Standard deviation of random perturbation 

% Seismicity catalogue 'cat'
cat = get_seismicity_sample_from_Dieterich94( ...
    neq, ...
    n0, e0, d0, ...
    stk, dip, ...
    cradius, ...
    dx); 

hf4 = figure(204); clf; hold on; grid on; box on; axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
%lim = get_cubic_domain_XYZ(n0, e0, d0, ceil(max(cat.ri)) );
lim  = get_cubic_domain_XYZ(n0, e0, d0, 1.5*finsrc.length);
set(gca,'xlim', lim.x, ...
        'ylim', lim.y, ...
        'zlim', lim.z, ...
        'yDir','reverse', ...
        'zDir','reverse', ...
        'view', view_angle)
    

plot_rectangular_slippatch(finsrc, n0, e0, d0, 'm')
plot3(cat.n, cat.e, cat.d, '.k')

1;





%% Infer its orientation with Principal Component Analysis
EVD = get_eigenval_decomp_NED(cat);

plot_eigenvectors(EVD.V, 10)







%% Projection of FM vectors onto plane
% We may want to project e.g. the slip vector onto a plane, e.g. to make a
% cross-section through a seismogenic volume

% To project a *POINT* onto a plane, we can use the following function
P1     = [0,  7, 0]; 
P0     = [0, -2, 0]; 
strike = 40;        
dip    = 70;         
plotme = true;

[P1p, r, rS, rD] = project_point_onto_plane(P1, P0, strike, dip, plotme);


% To project a *VECTOR* onto a plane, we can do the same thing, for the
% starting and end points of the vector. E.g. for the slip vector: 

s1_beg = [n0, e0, d0];         % Slip vector start point
s1_end = [n0+fm.s1(1), ...    % Slip vector end   point
          e0+fm.s1(2), ...
          d0+fm.s1(3)];

[s1_beg_p, r, rS, rD] = project_point_onto_plane(s1_beg, P0, strike, dip);
[s1_end_p, r, rS, rD] = project_point_onto_plane(s1_end, P0, strike, dip);
      


% Create a synthetic seismicity cloud
neq = 1e2;
xx  = 10 + 1.0*randn(neq,1);
yy  = 10 + 1.1*randn(neq,1);
zz  = 10 + 2.0*randn(neq,1);
points = [xx, yy, zz];

[P_projected, r, rS, rD] = project_pointcloud_onto_plane(points, ...
    P0, ...
    strike, ...
    dip, ...
    plotme);


% Add FMs to plot slip vectors

























%% APPENDIX
% Functions below are interesting, but not yet checked
if false
    
    %% Second plane
    [stk2,dip2,rak2,~,~] = focal_pl2pl(stk, dip, rak); % Compute second FM plane
    %[stk3,dip3,rak3]     = sdr2sdr    (stk, dip, rak); % NDs solution
    plot_beachball_pair(stk,dip,rak,stk2,dip2,rak2)
    % beachball(stk2,dip2,rak2)
    % plot_n_beachballs(stk,dip,rak,stk2,dip2,-90)
    
    
    
    
    
    %% FM angles to slip vector, and back
    [n1,n2,n3,s1,s2,s3]          = fm2nsb_vectors(stk,dip,rak);    % Turn FM angles into slip- and normal vectors
    [stkX,dipX,rakX,dipdir,ierr] = focal_nd2pl(n1,n2,n3,s1,s2,s3); % Turn slip- and normal vectors inot FM angles
    [rotangle,~,~]               = kagan([stk,dip,rak],[stkX,dipX,rakX]);
    
    % Given strike, dip and rake of a fault plane calculate the slip vector.
    [paz,pdip,taz,tdip,baz,bdip] = sdr2ptb(stk,dip,rak);
    
    
    
    %% Rotate FM
    
    % Randomly sample rotation angle from Von Mises distribution
    npu      = 30; % Nodal plane uncertainty
    oneSigma = npu*pi/180;
    kappa    = 1/oneSigma^2;             % kappa =! 1/sigma^2 [1/rad]
    rota0    = 90*pi/180;
    rota     = circ_vmrnd(rota0, kappa, 1e4)*180/pi;
    clf; hold on; histogram(rota)
    prctile(rota,[16,84])
    
    
    % Rotate FM around randomly chosen axis
    rota = 14;
    [stk2,dip2,rak2] = FMrot_known_axis (stk,dip,rak,rota,0,0,1);
    [stk3,dip3,rak3] = FMrot_random_axis(stk,dip,rak,rota);
    [stk,dip,rak; stk2,dip2,rak2; stk3,dip3,rak3];
    
    % Check if angle is right
    [rotangle,~,~] = kagan([stk2,dip2,rak2],[stk,dip,rak]);
    rotangle-rota
    
    
    
    %% Measure FM rotate angle
    %  1. Considering FP ambiguity
    [rota1,~,~] = kagan([stk,dip,rak],[stk2,dip2,rak2]);
    
    %  2. NOT considering FP ambiguity (e.g. for finding preferential fault
    %  plane)
    rota2 = fault_plane_rot(stk,dip,rak,stk2,dip2,rak2);
    
    
    study_rotation_angle_and_FP_ambiguity
    
    
    
    
    % Rotation angle between main and auxiliary planes
    % with    ambiguity: 0°
    % without ambiguity: 180°
    stk = 80;
    dip = 20;
    rak = 40;
    
    [stk2,dip2,rak2,~,~] = focal_pl2pl(stk, dip, rak); % Compute second FM plane
    [rota1,~,~]          = kagan([stk,dip,rak],[stk2,dip2,rak2]);
    rota2                = fault_plane_rot(stk,dip,rak,stk2,dip2,rak2);
    
    
    
    
    
    %% Fault plane ambiguity
    study_fault_plane_ambiguity
    
    
    
    
    %% Find preferred FP based on similarity with reference mechanism
    %  Compute 2nd fault plane and choose the one that is closer to ref FM
    %  Not sure if this includes FP ambiguity or not
    [str,dip,rak] = find_pref_faultMech(str,  dip,  rak, ...
        str0, dip0, rak0);
    
    
    
    
    
    
    %% Slip vectors
    
    % For normal FM with strike = 0, slip vector points to right, i.e. has
    % n2>0; for same mechanism with strike = 180, sign of n2 is flipped
    stk =   0;
    dip =  45;
    rak = -90;
    % beachball(stk, dip, rak)
    
    % Turn FM params into slip- and normal vectors [deg]
    [n1,n2,n3,s1,s2,s3] = fm2nsb_vectors(stk,dip,rak);
    
    [stk2,dip2,rak2,~,~] = focal_pl2pl(stk, dip, rak); % Compute second FM plane
    [n1,n2,n3,s1,s2,s3]  = fm2nsb_vectors(stk2,dip2,rak2);
    cat
    
    
    
end

%% APPENDIX
%[stk2,dip2,rak2] = ComputeSecondPlane(stk,dip,rak);