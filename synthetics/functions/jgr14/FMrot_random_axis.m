function [stk,dip,rak] = FMrot_random_axis(stk0,dip0,rak0,rota) 
% Rotate FM around randomly oriented rotation axis
% UNITS: all input-angles are to be given in [DEG]

% Find a random unit vector with the method of Marsaglia (1972) 
[r1,r2,r3] = sphere_point_picking_Marsaglia;

% Turn FM params into slip- and normal vectors
% Input-units: [deg]
[n1,n2,n3,s1,s2,s3] = fm2nsb_vectors(stk0,dip0,rak0);

% Rotate both, slip- and normal vectors around the rotation axis r by <rota>
% degrees, corresponding to the catalog FM error
[n1r,n2r,n3r] = rodrigues_rotation(n1,n2,n3,r1,r2,r3,rota);
[s1r,s2r,s3r] = rodrigues_rotation(s1,s2,s3,r1,r2,r3,rota);

% Turn rotated slip- and normal vectors back to FMs
% Input-units: [deg]
% Output-units: [deg]
[stk,dip,rak,dipdir,ierr] = focal_nd2pl(n1r,n2r,n3r,s1r,s2r,s3r);


% stk = round(stk);
% dip = round(dip);
% rak = round(rak);
if stk==360; stk=0; 
end