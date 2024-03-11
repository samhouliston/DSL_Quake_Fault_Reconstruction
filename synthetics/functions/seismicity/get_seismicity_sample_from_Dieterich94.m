function cat = get_seismicity_sample_from_Dieterich94( ...
    neq, ...
    n0, e0, d0, ...
    azi, dip, ...
    crack_radius, ...
    dx)


% Get samples of quakes from stress profile
[cat, stress] = sample_Dieterich94_stress_profile(...
    neq, ...
    crack_radius, ...
    1, ...
    0);


% Distribute samples along plane with particular dip and rake
% Use x- and y-values from sampled seismicity as along-strike and along-dip
% distances, respectively
[n, e, d] = azidip2plane_NED(n0, e0, d0, ...
                             azi, dip, ...
                             cat.xi, ...
                             cat.yi);

% clf; hold on; axis equal; % plot3(e, n, h, 'ok')

% Generated samples are on plane. Add random perturbation to simulate 
% off-plane seismicity
cat.n = n +dx*randn(neq,1);
cat.e = e +dx*randn(neq,1);
cat.d = d +dx*randn(neq,1);

cat.azi    = azi;
cat.dip    = dip;
cat.stress = stress;

cat.prop.tstring = sprintf('Seismicity sampled from Dieterich 1994 stress profile, for crack with radius = %im', crack_radius);                     
cat.prop.name    = 'Dieterich 1994 seismicity';









 








%% Appendix


%as   = get_seismicity_sample_from_Dieterich94(nrnd,crad)

%[as,stress] = sample_Dieterich94xm1_stress_profile(nrnd,crad); % Dieterich 1994 but with r^-1 stress decay

% rp.e = rp.e0 +dxrand*randn(size(rp.e));
% rp.n = rp.n0 +dxrand*randn(size(rp.n));
% rp.h = rp.h0 +dxrand*randn(size(rp.h));
% plot3(rp.e, rp.n, rp.h,'o', ...
%     'markerSize',5, ...
%     'markerFaceColor',[.9 0 0],...
%     'markerEdgeColor',[.7 .7 .7])


% dtau0       = 10;
% dtaumax     = 10;
% [as,stress] = sample_centerxm1_stress_profile(nrnd,dtau0,dtaumax); % r^-1 stress decay from center of crack
%                                                                    % Doesn't have a crack-length, so can't be evaluated...
% 
% % % r^-1 stress decay from center of crack; non-circular and rotated
% nrnd        = 4e2;
% dtau0       = 10;
% dtaumax     = 10;
% [as,stress] = sample_rotatedEllipsexm1_stress_profile(nrnd,dtau0,dtaumax);
% as.dx = as.dx*3;
% as.dy = as.dy +.4*as.dx;
% clf; plot(as.dx,as.dy,'.k')
% axis equal