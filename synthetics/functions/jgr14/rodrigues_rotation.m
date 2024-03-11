function [v1_rot,v2_rot,v3_rot] = rodrigues_rotation(v1,v2,v3,r1,r2,r3,rota)
% rotates an input vector V=(v1,v2,v3) about a rotation axis R=(r1,r2,r3) by 
% UNITS: rota in [DEG]
% MAM, 110207

deg2rad = pi()/180;
theta = rota*deg2rad;

% Scale rotation axis to unit length
scale = sqrt(r1^2+r2^2+r3^2);
r1_u    = r1/scale;
r2_u    = r2/scale;
r3_u    = r3/scale;
%disp(num2str(sqrt(r1_u^2+r2_u^2+r3_u^2)))

% % Crossproduct of V and R
% VcR_1 = v2*r3_u - v3*r2_u;
% VcR_2 = v3*r1_u - v1*r3_u;
% VcR_3 = v1*r2_u - v2*r1_u;

% Crossproduct of R and V --> the one to use!
RcV_1 = r2_u*v3 - r3_u*v2;
RcV_2 = r3_u*v1 - r1_u*v3;
RcV_3 = r1_u*v2 - r2_u*v1;

% Rotate V
v1_rot = v1*cos(theta)+RcV_1*sin(theta)+r1_u*(v1*r1_u+v2*r2_u+v3*r3_u)*(1-cos(theta));
v2_rot = v2*cos(theta)+RcV_2*sin(theta)+r2_u*(v1*r1_u+v2*r2_u+v3*r3_u)*(1-cos(theta));
v3_rot = v3*cos(theta)+RcV_3*sin(theta)+r3_u*(v1*r1_u+v2*r2_u+v3*r3_u)*(1-cos(theta));

% If one of the components is close to 0, set it to zero
if ((v1_rot < 1E-5) && (v1_rot > -1E-5))
    v1_rot = 0;end
if ((v2_rot < 1E-5) && (v2_rot > -1E-5))
    v2_rot = 0;end
if ((v3_rot < 1e-05) && (v3_rot > -1E-5))
    v3_rot = 0;end