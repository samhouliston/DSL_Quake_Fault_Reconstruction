function catFM = get_finite_sources_from_FM_cat(catFM, stressdrop)

nfm = numel(catFM.mag);
catFM.finsrc = cell(nfm,1);
for ifm = 1:nfm
    
    print_iter_nums(ifm,nfm,1000)
    
    src = get_finite_source_from_FM(catFM.lat(ifm), catFM.lon(ifm), catFM.dep(ifm), ...
                                    catFM.stk(ifm), catFM.dip(ifm), catFM.mag(ifm), ...
                                    stressdrop);
	catFM.finsrc{ifm} = src;
end