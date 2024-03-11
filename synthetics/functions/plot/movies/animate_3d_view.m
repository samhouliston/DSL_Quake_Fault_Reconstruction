function animate_3d_view(az,dip, videoFullName)
% n   = 100;
% az  = linspace(320, 680, n)';
% dip = linspace( 50, 5,n)';

if nargin<3; videoFullName = [];
             printme       = false;
else         printme       = true;
end

if printme 
    V = VideoWriter(videoFullName, 'MPEG-4');
    V.FrameRate=8; % Frames per second
    open(V);
end


nn = numel(az);
for in = 1:nn
    
    set(gca,'view',[az(in) dip(in)])
    pause(.1)
    
    %figure('units','pixels','position',[0 0 1440 1080])
    
    if printme
        frame = getframe(gcf);
        writeVideo(V, frame);
    end
end

if printme
    close(V)
end
