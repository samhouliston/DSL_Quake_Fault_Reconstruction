cellFault =   {[-116.449, 34.362;
                -116.421, 34.094;
                -116.421, 34.094;
                -116.449, 34.362;
                -116.449, 34.362],...
               [-116.545, 34.524;
                -116.416, 34.306;
                -116.416, 34.306;
                -116.545, 34.524;
                -116.545, 34.524],...
               [-116.716, 34.696;
                -116.464, 34.448;
                -116.464, 34.448;
                -116.716, 34.696;
                -116.716, 34.696]};

SHP_FILE    = 'CFM5.2_shp\CFM52_noZ_noB.shp';
bbox        = [-116.75 34 ;-116.4 35];
sel     = [1 9 10 13:19 25 27];
%
shp     = shaperead(SHP_FILE,'BoundingBox',bbox, 'UseGeoCoords', true);
shp     = shp(:);

SHP_FILE    = 'CFM5.2_shp\CFM52_noZ.shp';
sel         = [3:4];
shp2     = shaperead(SHP_FILE,'BoundingBox',bbox, 'UseGeoCoords', true);
shp2     = shp2(:);

shpAll      = [shp; shp2];
cellLtLn    = cell(numel(shpAll),1);
figure;
for i=1:numel(shpAll)
    [x,y] = captopUTM(shpAll(i).Lat(1:end-1)',shpAll(i).Lon(1:end-1)');
    cellLtLn{i} = [x-mX y-mY]/1000;
    %cellLtLn{i} = [x y];
    plot(cellLtLn{i}(:,1), cellLtLn{i}(:,2));
    hold on;
end



geoshow(shp);
for i=1:numel(shp)
    figure;
    plot(shp(i).Lon,shp(i).Lat);
    title([num2str(i) ': ' shp(i).Layer])
    xlim([-117 -116.25]);
    ylim([33.8 34.8]);
    daspect([1 1 1]);
end