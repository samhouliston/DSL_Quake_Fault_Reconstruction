function lim = get_bounding_box(XYZ, plo, pup)

% % Include all points between [plo pup]
% plo = 50-p/2; % Upper percentile
% pup = 50+p/2; % Lower percentile

if nargin<2; 
    plo = 0;
    pup = 100;
end

x     = XYZ(:,1);
y     = XYZ(:,2);
lim.x = [prctile(x, plo), prctile(x, pup)];
lim.y = [prctile(y, plo), prctile(y, pup)];

lim.x = sort(lim.x, 'ascend');
lim.y = sort(lim.y, 'ascend');

ndim = size(XYZ,2);
if ndim==3
    z     = XYZ(:,3);
    lim.z = [prctile(z, plo), prctile(z, pup)];
    lim.z = sort(lim.z, 'ascend');
end