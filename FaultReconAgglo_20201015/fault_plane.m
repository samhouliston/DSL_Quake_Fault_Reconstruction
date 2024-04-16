function rot_pnt = fault_plane(ST,DP,L,W,D,dens,rseed,loc_err,MODE_ELIP, DO_PLOT)
% Sample a single fault plane with random points
% ST        : Strike-angle from North (in degrees)
% DP        : Dip-angle (in degrees, must be scalar)
% L         : Fault length in strike direction (LEN > 0)
% W         : Fault width in dip direction (WIDTH > 0)
% D         : Depth of the fault centroid (DEP > 0)
% dens      : Event density per km2
% rseed     : Random seed generator for repeatable synthetics
% loc_err   : Thickness of the fault (along the smallest eigenvector)
% MODE_ELIP : 0: Rectangular planes, 1: Gaussian ellipses 
% DO_PLOT   : 0: no plot, 1: plot

% Y.Kamer 20201015

rng(rseed); 
npo     = round(dens*L*W);
%loc_err = 1; %thickness * 2
strike  = ST*pi/180;
dip     = DP*pi/180;
strike2  = (90-ST)*pi/180;
dip2    = (DP-90)*pi/180;
depth   = D;
d       = depth + sin(dip).*W/2;

    alpha = pi/2 - strike;
	x_fault = L/2*cos(alpha)*[-1,1,1,-1] + sin(alpha)*cos(dip)*W/2*[-1,-1,1,1];
	y_fault = L/2*sin(alpha)*[-1,1,1,-1] + cos(alpha)*cos(dip)*W/2*[1,1,-1,-1];
	z_fault = -d + sin(dip)*W*[1,1,0,0];

    x_pnt1 = (rand(npo,1)-0.5)*L;
    y_pnt1 = (rand(npo,1)-0.5)*0.1;
    z_pnt1 = (rand(npo,1)-0.5)*W;
    
    
 
    x_pnt = L*cos(alpha)*x_pnt1 + sin(alpha)*cos(dip)*W*y_pnt1;
    y_pnt = L*sin(alpha)*x_pnt1 + cos(alpha)*cos(dip)*W*y_pnt1;
    z_pnt = sin(dip)*W*z_pnt1;
    
    %{
    t=dip;
    dip=alpha;
    alpha=t;
    %}
    v     = [0 0 1];
    mat_r = rotMat(v,strike2);
    
    rot_pnt = [x_pnt1  y_pnt1 z_pnt1]*mat_r';
    rot_v   = [1 0 0]*mat_r';
    
    mat_r = rotMat(rot_v,dip2);
    rot_pnt = rot_pnt*mat_r';
    rot_pnt(:,3)=rot_pnt(:,3)-depth;
    
    rot_pnt = rot_pnt+(randn(size(rot_pnt)))*loc_err/2;
    if(MODE_ELIP) %fit the planes a gaussian ellipsoid and sample the points from it
        rot_pnt = mvnrnd(mean(rot_pnt),cov(rot_pnt),size(rot_pnt,1));
    end
    if(DO_PLOT)
        %patch(x_fault,y_fault,z_fault,.3*[1,1,1],'EdgeColor','k','LineWidth',2)
        hold on;
        %scatter3(x_pnt,y_pnt,z_pnt);
        scatter3(rot_pnt(:,1),rot_pnt(:,2),rot_pnt(:,3));
        view(3);
        grid on;
        xlabel('X');ylabel('Y');zlabel('Z');
        daspect([1 1 1]);
    end
    
    function [mat_R]=rotMat(v,theta)
        v = v'; 
        v_mat = [0 -v(3) v(2); v(3) 0 -v(1);-v(2) v(1) 0];
        mat_R = eye(3)*cos(theta)+sin(theta)*v_mat+(1-cos(theta))*(v*v');
    end
end

