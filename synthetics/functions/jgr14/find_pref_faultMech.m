function [strike,dip,rake]=find_pref_faultMech(strike1,dip1,rake1, str_master,dip_master,rak_master)
% Find the one fault mechanism (= uniquely defined slip and normal vectors) 
% of a FM which is closer to a "master" fault mechanism
%
% Units: all in- and output-angles are in [DEG]

[strike2,dip2,rake2] = ComputeSecondPlane(strike1,dip1,rake1);

strike2 = round(strike2);
dip2    = round(dip2);
rake2   = round(rake2);

if strike2==360; strike2=0; end

% Determine which nodal plane is closer to the
% regional master mechanism by computing the minimal rotation angles.
[rota1] = fault_plane_rot(strike1,dip1,rake1,str_master,dip_master,rak_master);
[rota2] = fault_plane_rot(strike2,dip2,rake2,str_master,dip_master,rak_master);

% addpath(genpath('~/programs/matlab/nDeichmann'))
% clf; axis equal
% beachball_mod(strike1,dip1,rake1,1,1,0.3,1)
% beachball_mod(strike2,dip2,rake2,2,1,0.3,1)

% If FM2 is closer to the regional master FM, use it instead of the catalog FM
if     rota2>=rota1; strike = strike1;
                     dip    = dip1;
                     rake   = rake1;
elseif rota2< rota1; strike = strike2;
                     dip    = dip2;
                     rake   = rake2;
end