function param_All=paramJoin(param_All,param_Inc)
    numALL              = numel(param_All.w);
    numINC              = numel(param_Inc.w);
    param_All.m         = [param_All.m param_Inc.m];
    param_All.w         = [param_All.w param_Inc.w];
    param_All.bkg       = [param_All.bkg param_Inc.bkg];
    param_All.krn_type  = [param_All.krn_type param_Inc.krn_type];
    param_All.sphR      = [param_All.sphR param_Inc.sphR];
    param_All.sig       = [param_All.sig param_Inc.sig];
    for i=1:numINC
        param_All.covar(:,:,numALL+i)   = param_Inc.covar(:,:,i);
        param_All.bbox(:,:,numALL+i)    = param_Inc.bbox(:,:,i);
    end
    %{
    IDbkg   = find(param_All.bkg);
    if(numel(IDbkg)==2)
         vert1 = DrawCuboidYK(param_All.m(:,IDbkg(1)),param_All.covar(:,:,IDbkg(1)));
         vert2 = DrawCuboidYK(param_All.m(:,IDbkg(2)),param_All.covar(:,:,IDbkg(2)));
         vertJ = [vert1 vert2];
         [boxCen, covMat]           = minboundboxYK(vertJ(1,:),vertJ(2,:),vertJ(3,:),'v',1);
         vol1   = prod(sqrt(eig(param_All.covar(:,:,IDbkg(1)))));
         vol2   = prod(sqrt(eig(param_All.covar(:,:,IDbkg(2)))));
         volJ   = vol1+vol2;
         volBnd = prod(sqrt(eig(covMat)));
         volMrg = min(volJ,volBnd);
         sclVol = (volMrg/volBnd)^(2/3);
         param_All.m(:,IDbkg(1))    = boxCen;
         param_All.w(IDbkg(1))      = param_All.w(IDbkg(1))+param_All.w(IDbkg(2));
         param_All.sphR(IDbkg(1))   = param_All.sphR(IDbkg(1))+param_All.sphR(IDbkg(2));
         param_All.covar(:,:,IDbkg(1))=covMat*sclVol;
         
        param_All.m(:,IDbkg(2))         = [];
        param_All.covar(:,:,IDbkg(2))   = [];
        param_All.w(IDbkg(2))           = [];
        param_All.bkg(IDbkg(2))         = [];
        param_All.sphR(IDbkg(2))        = [];
        param_All.sig(:,IDbkg(2))       = [];
        param_All.krn_type(IDbkg(2))    = [];
    end
    %}
end