function [Strike_rot,Dip_rot,Rake_rot] = FMrot(Strike,Dip,Rake,rota) 
% Rotate FM around randomly oriented rotation axis
% UNITS: all input-angles are to be given in [DEG]

% Find a random unit vector with the method of Marsaglia (1972) 
[r1,r2,r3] = sphere_point_picking_Marsaglia;

% Turn FM params into slip- and normal vectors
% Input-units: [deg]
[n1,n2,n3,s1,s2,s3] = fm2nsb_vectors(Strike,Dip,Rake);

% Rotate both, slip- and normal vectors around the rotation axis r by <rota>
% degrees, corresponding to the catalog FM error
[n1_rot,n2_rot,n3_rot] = rodrigues_rotation(n1,n2,n3,r1,r2,r3,rota);
[s1_rot,s2_rot,s3_rot] = rodrigues_rotation(s1,s2,s3,r1,r2,r3,rota);

% Turn rotated slip- and normal vectors back to FMs
% Input-units: [deg]
% Output-units: [deg]
[Strike_rot,Dip_rot,Rake_rot,dipdir,ierr] = focal_nd2pl(n1_rot,n2_rot,n3_rot,s1_rot,s2_rot,s3_rot);
