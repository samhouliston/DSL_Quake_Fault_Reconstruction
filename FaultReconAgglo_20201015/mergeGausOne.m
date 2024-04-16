function newPrm=mergeGausOne(prm,vi,vj,KEEP_WGT)
%if KEEP_WGT is provided keep the calculated merged weight, instead of 1
    for k=1:numel(vi)
        i       = vi(k);
        j       = vj(k);
        w       = prm.w(i) + prm.w(j);
        if(nargin>3)
            newPrm.w(k)     = w;
        else
            newPrm.w(k)     = 1;
        end
        %
        if(~prm.bkg(i) && ~prm.bkg(j))
            BKG_NUM = 0; %none is bkg
        elseif(prm.bkg(i) && prm.bkg(j))
            BKG_NUM = 2; %both are bkg
            jSC     = 1;
            iSC     = 1;
        elseif(prm.bkg(i))
            BKG_NUM = 1; %only i is bkg
            jSC     = 8;
            iSC     = 1;
        elseif(prm.bkg(j))
            BKG_NUM = 1; %only j is bkg
            jSC     = 1;
            iSC     = 8;
        end
        
        newPrm.mix_type             = 'gaus';
        newPrm.krn_type(:,k)        = 2;
        newPrm.sphR(k)              = max(prm.sphR(i),1)+max(prm.sphR(j),1);
        vertJ                       = [prm.bbox(:,:,i);prm.bbox(:,:,j)]';
        [corPts, boxCen, covMat]    = minboundboxYK(vertJ(1,:),vertJ(2,:),vertJ(3,:),'v',1);
        newPrm.bbox(:,:,k)          = corPts;       
        switch BKG_NUM
            case 0
                newPrm.bkg(k)         = false;
                m       = (1/w)*(prm.w(i)*prm.m(:,i) + prm.w(j)*prm.m(:,j));
                covar   = (prm.w(i)/w)*(prm.covar(:,:,i) + (prm.m(:,i)-m)*(prm.m(:,i)-m)')...
                        + (prm.w(j)/w)*(prm.covar(:,:,j) + (prm.m(:,j)-m)*(prm.m(:,j)-m)');
                newPrm.m(:,k)         = m;
                newPrm.sig(:,k)       = sqrt(diag(covar));
                newPrm.covar(:,:,k)   = covar;
                %newPrm.bbox(:,:,k)    = DrawCuboidYK(m,covar*12)';  
            case {1,2}
                newPrm.bkg(k)       = true;
                %{
                vert1 = DrawCuboidYK(prm.m(:,i),prm.covar(:,:,i)*iSC);
                vert2 = DrawCuboidYK(prm.m(:,j),prm.covar(:,:,j)*jSC);
                vertJ = [vert1 vert2];
                [boxCen, covMat] = minboundboxYK(vertJ(1,:),vertJ(2,:),vertJ(3,:),'v',1);
                %}
                %{
                volBnd = prod(sqrt(eig(covMat)));
                vol1   = prod(sqrt(eig(prm.covar(:,:,i))));
                vol2   = prod(sqrt(eig(prm.covar(:,:,j))));
                volJ   = vol1+vol2;
                volMrg = min(volJ,volBnd);
                sclVol = (volMrg/volBnd)^(2/3);
                newPrm.sig(:,k)       = sqrt(diag(covMat*sclVol));
                newPrm.covar(:,:,k)   = covMat*sclVol;
                %}
                newPrm.m(:,k)         = boxCen;
                newPrm.sig(:,k)       = sqrt(diag(covMat));
                newPrm.covar(:,:,k)   = covMat;
        end
        %}
        %{
        if(prm.bkg(i))
            IDbkg   = i;
            bkg     = true;
        elseif(prm.bkg(j))
            IDbkg   = j;
            bkg     = true;
        else
            bkg     = false;
        end 
        newPrm.bkg(k)         = bkg;
        
        if(bkg)
            newPrm.m(:,k)         = prm.m(:,IDbkg);
            newPrm.sig(:,k)       = prm.sig(:,IDbkg);
            newPrm.covar(:,:,k)   = prm.covar(:,:,IDbkg);
        else
        
            m       = (1/w)*(prm.w(i)*prm.m(:,i) + prm.w(j)*prm.m(:,j));
            covar   = (prm.w(i)/w)*(prm.covar(:,:,i) + (prm.m(:,i)-m)*(prm.m(:,i)-m)')...
                    + (prm.w(j)/w)*(prm.covar(:,:,j) + (prm.m(:,j)-m)*(prm.m(:,j)-m)');
            newPrm.m(:,k)         = m;
            newPrm.sig(:,k)       = sqrt(diag(covar));
            newPrm.covar(:,:,k)   = covar;
        end
        %}

    end
end