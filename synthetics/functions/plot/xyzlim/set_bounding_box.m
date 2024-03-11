function lim = set_bounding_box(XYZ, plo, pup)

% Find and set x/y/z limits to include all data points that lie between
% the <plo> and <pup> pecentiles

% Examples: include all data points from 1st to 99th percentile
% 2D: set_bounding_box([x.val,y.val], 1, 99)
% 3D: set_bounding_box([x.val,y.val,z.val], 1, 99)

if nargin<2; 
    plo = 0;
    pup = 100;
end

ndim = size(XYZ,2);

lim = get_bounding_box(XYZ, plo, pup);

if ndim==2
	set(gca, 'xlim', lim.x, ...
             'ylim', lim.y)

elseif ndim==3
    set(gca, 'xlim', lim.x, ...
             'ylim', lim.y, ...
             'zlim', lim.z) 

else
    error('Works only in 2D or 3D')
end

lim.plo = plo;
lim.pup = pup;