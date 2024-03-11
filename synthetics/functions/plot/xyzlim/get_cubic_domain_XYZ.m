function lim = get_cubic_domain_XYZ(x0, y0, z0, sidelength)
% compute the x/y/z limits for a cube with its centre at x0/y0/z0

dx = sidelength/2;

lim.x = [x0-dx x0+dx];
lim.y = [y0-dx y0+dx];
lim.z = [z0-dx z0+dx];