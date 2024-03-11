function phi = get_fault_plane_rot(stk1,dip1,rak1,stk2,dip2,rak2)

% Fault plane rotation: computes the minimum rotation angle (rota) between 
% two slip mechanisms, i.e. between two sets of normal and slip vectors,
% without considering any fault plane ambiguity.
%
% Sign conventions are from Aki and Richards
% All in- and outputs in DEGREES
%
% I *think* this is from an orignal script provided by JHardebeck in 2013,
% Modified by Men-Andrin Meier in the same year, and now with minor
% adaptions in 2022
%
% Some interesting examples that illustrate the difference between this
% simple rotation and a conventional FM rotation that heeds the FM
% ambiguity:
% [rota]=fault_plane_rot(90,90,-180,0,90,0);
% [rota]=fault_plane_rot(90,90,-180,0,90,-180);

deg2rad = pi/180;
rad2deg = 180/pi;



%% fm2tpb - turn FMs into T- P- and B-axes
% Formulae taken from Kagan 2007

% FM1 --------------------------------------------------
stk1 = stk1*deg2rad; 
dip1 = dip1*deg2rad; 
rak1 = rak1*deg2rad;

t1_1 = (-sin(stk1)*sin(dip1) + cos(stk1)*cos(rak1) + sin(stk1)*cos(dip1)*sin(rak1))/sqrt(2);
t1_2 = (cos(stk1)*sin(dip1) + sin(stk1)*cos(rak1) - cos(stk1)*cos(dip1)*sin(rak1))/sqrt(2);
t1_3 = (-cos(dip1) - sin(dip1)*sin(rak1))/sqrt(2);
t1   = [t1_1; t1_2; t1_3];

p1_1 = (-sin(stk1)*sin(dip1) - cos(stk1)*cos(rak1) - sin(stk1)*cos(dip1)*sin(rak1))/sqrt(2);
p1_2 = (cos(stk1)*sin(dip1) - sin(stk1)*cos(rak1) + cos(stk1)*cos(dip1)*sin(rak1))/sqrt(2);
p1_3 = (-cos(dip1) + sin(dip1)*sin(rak1))/sqrt(2);
p1   = [p1_1; p1_2; p1_3];

b1_1 = cos(stk1)*sin(rak1) - sin(stk1)*cos(dip1)*cos(rak1); 
b1_2 = sin(stk1)*sin(rak1) + cos(stk1)*cos(dip1)*cos(rak1); 
b1_3 = sin(dip1)*cos(rak1);
b1   = [b1_1; b1_2; b1_3];

% Check for non-orthogonal b-vector
%b1_orth = cross(t1,p1);
%db=b1-b1_orth

TPB1     = [t1, p1, b1];
TPB1_inv = [-t1, -p1, b1];




% FM2  --------------------------------------------------
stk2 = stk2*deg2rad; 
dip2 = dip2*deg2rad; 
rak2 = rak2*deg2rad;

t2_1 = (-sin(stk2)*sin(dip2) + cos(stk2)*cos(rak2) + sin(stk2)*cos(dip2)*sin(rak2))/sqrt(2);
t2_2 = ( cos(stk2)*sin(dip2) + sin(stk2)*cos(rak2) - cos(stk2)*cos(dip2)*sin(rak2))/sqrt(2);
t2_3 = (-cos(dip2) - sin(dip2)*sin(rak2))/sqrt(2);
t2   = [t2_1; t2_2; t2_3];

p2_1 = (-sin(stk2)*sin(dip2) - cos(stk2)*cos(rak2) - sin(stk2)*cos(dip2)*sin(rak2))/sqrt(2);
p2_2 = ( cos(stk2)*sin(dip2) - sin(stk2)*cos(rak2) + cos(stk2)*cos(dip2)*sin(rak2))/sqrt(2);
p2_3 = (-cos(dip2) + sin(dip2)*sin(rak2))/sqrt(2);
p2   = [p2_1; p2_2; p2_3];

b2_1 = cos(stk2)*sin(rak2) - sin(stk2)*cos(dip2)*cos(rak2); 
b2_2 = sin(stk2)*sin(rak2) + cos(stk2)*cos(dip2)*cos(rak2); 
b2_3 = sin(dip2)*cos(rak2);
b2   = [b2_1; b2_2; b2_3];

TPB2 = [t2, p2, b2];
TPB2_inv = [-t2, -p2, b2];

% Some (non-systematic) checks for orthogonality
test1=t1(3)^2+p1(3)^2+b1(3)^2;
test2=t2(3)^2+p2(3)^2+b2(3)^2;
test3=t1(1)^2+t1(2)^2+t1(3)^2;
test=(test1+test2+test3)/3;
%fprintf(1,['orthogonality check: ', num2str(test), '\n'])
% has to be equalt to 1

%% Compute rotation matrix R between the 2 matrices and rot-angle phi
R     = TPB1*TPB2';
R_inv = TPB1_inv*TPB2';
% There are two more possible combinations: 
% TPB2*TPB1_inv' and TPB1_inv*TPB2_inv' but they are equivalent to the
% above two.

% Determine rotation angle phi, the smaller of which is the minimum
% rotation anlge between the two slip mechanisms.
phi1 = acos((R(1,1)+R(2,2)+R(3,3)-1)/2);
phi1 = phi1*rad2deg;

phi2 = acos((R_inv(1,1)+R_inv(2,2)+R_inv(3,3)-1)/2);
phi2 = phi2*rad2deg;

if (phi1 < phi2)
    phi = phi1;
else
    phi = phi2;
end