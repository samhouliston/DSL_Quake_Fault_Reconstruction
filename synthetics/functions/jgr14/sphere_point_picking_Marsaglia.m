function [x,y,z] = sphere_point_picking_Marsaglia
% Finds a truely random point on a sphere (= a truely random unit vector), 
% such that if executed a large number of times, the points will be 
% uniformely distributed over the sphere surface.
% MAM 110208

repeat = 1;

while repeat
x1=-1+2*rand(1,1);
x2=-1+2*rand(1,1);

if ((x1*x1 + x2*x2) < 1)
    repeat = 0;
end

end

% The components of the random unit vector are then
x = 2*x1*sqrt(1-x1*x1-x2*x2);
y = 2*x2*sqrt(1-x1*x1-x2*x2);
z = 1-2*(x1*x1+x2*x2);