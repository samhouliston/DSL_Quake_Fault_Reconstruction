function plot_rectangular_slippatch(finsrc, n0, e0, d0, colour)

if nargin<5; colour = [0 0 1];
end

       
plot3(n0, e0, d0, 'ok', 'markerFaceColor', 'r')

fill3( finsrc.x5(:), finsrc.y5(:), finsrc.z5(:), ...
    ones(5,1), ...
    'faceColor', colour, ...
    'edgeColor', [.2 .2 .2], ...
    'faceAlpha', .3)

% lim = get_cubic_domain_XYZ(n0, e0, d0, 1.3*finsrc.length);
% set(gca,'xlim', lim.x, ...
%         'ylim', lim.y, ...
%         'zlim', lim.z)