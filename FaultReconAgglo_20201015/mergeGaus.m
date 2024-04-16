function prm=mergeGaus(prm,vi,vj)
    for k=1:numel(vi)
        i           = vi(k);
        j           = vj(k);
        numK        = numel(prm.w);
        newPrm      = mergeGausOne(prm,i,j,1);
        
        prm.bkg(numK+1)         = newPrm.bkg;
        prm.w(numK+1)           = newPrm.w;
        prm.m(:,numK+1)         = newPrm.m;
        prm.sig(:,numK+1)       = newPrm.sig;
        prm.krn_type(numK+1)    = newPrm.krn_type;
        prm.sphR(numK+1)        = newPrm.sphR;
        prm.covar(:,:,numK+1)   = newPrm.covar;
        prm.bbox(:,:,numK+1)    = newPrm.bbox;
        
    end
    prm.bkg([vi vj])            = [];
    prm.w([vi vj])              = [];
    prm.m(:,[vi vj])            = [];
    prm.sig(:,[vi vj])          = [];
    prm.krn_type([vi vj])       = [];
    prm.sphR([vi vj])           = [];
    prm.covar(:,:,[vi vj])      = [];
    prm.bbox(:,:,[vi vj])       = [];
end