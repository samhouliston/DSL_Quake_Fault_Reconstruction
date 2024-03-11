function make_rotating_gif(gif_name)

% DOES NOT YET WORK. WHY IS THIS SO COMPLICATED. 2023.
gif_name = "testAnimated.gif"; % Specify the output file name
for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,gif_name,"gif","LoopCount",Inf,"DelayTime",1);
    else
        imwrite(A,map,gif_name,"gif","WriteMode","append","DelayTime",1);
    end
end