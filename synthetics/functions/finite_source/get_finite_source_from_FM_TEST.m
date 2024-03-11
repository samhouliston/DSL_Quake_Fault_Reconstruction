%function get_finite_source_from_FM_TEST

hbm = plot_Amatrice16_basemap(sta, F, catH5, finsrc);

lat = catH5.lat(catH5.ms.norcia.i);
lon = catH5.lon(catH5.ms.norcia.i);
dep = catH5.dep(catH5.ms.norcia.i);
mag = catH5.mag(catH5.ms.norcia.i);
stk = 150; 
dip = 55; 
rak = -90;

src = get_finite_source_from_FM(lat, lon, dep, stk, dip, mag);

fill3( src.lat5(:), src.lon5(:), src.dep5(:), ones(5,1), ...
    'faceColor', 'b', ...
    'edgeColor', 'k', ...
    'faceAlpha', .3)
