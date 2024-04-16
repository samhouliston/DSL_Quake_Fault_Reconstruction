function TF=touchBBOX(prm,i,j)
    TF = 0;
    if(prm.bkg(i) || prm.bkg(j))
        return;
    end
     for k=1:2
        if(k==1)
           a=i;b=j;
        else
           a=j;b=i;
        end
        verts           = prm.bbox(:,:,a);
        data            = prm.bbox(:,:,b);
        tol             = 1.e-13*mean(abs(verts(:)));
        inBKG           = inhull(data,verts,[],tol);
        TF  = TF || any(inBKG);
        if(TF)
            return;
        end
     end
    
end
