function iCond = param2iCond( param_best )
    SZ  = size(param_best.w);
    if(SZ(1)>SZ(2))
        param_best.w = param_best.w';
    end
    nK      = 1;
    iCond   = struct('centrCell',cell(nK,1),...
                    'covarCell',cell(nK,1),...
                    'bkgCell',cell(nK,1),...
                    'weightCell',cell(nK,1));
    iCond(1).centrCell{1}   = param_best.m';
    iCond(1).covarCell      = param_best.covar;
    iCond(1).weightCell     = param_best.w;
    if(isfield(param_best,'w'))
        iCond(1).bkgCell     = param_best.bkg;
    end
end

