function val_pdf = mvncdf_NOmonte(data,param,j,dataT)
    if(max(isnan(param.w)))
        val_pdf=ones(size(data,1),1)*0;
        return;
    end
    val_pdf=ones(size(data,1),1)*0;
    try
        switch param.mix_type
            case {'gaus','gausB'}
                
                if(param.bkg(j))
                    %verts           = DrawCuboidYK(param.m(:,j),param.covar(:,:,j));
                    %inBKG           = inhull(data,verts');
                    verts           = param.bbox(:,:,j);
                    tol             = 1.e-10*mean(abs(verts(:)));
                    inBKG           = inhull(data,verts,[],tol);
                    val_pdf(inBKG)  = param.w(j)/prod(sqrt(eig(param.covar(:,:,j))));

                else
                    if(nargin==4)
                        val_pdf = gaussianValue(dataT,param.m(:,j),param.covar(:,:,j))*param.w(j);
                    else
                        val_pdf = mvnpdf(data,param.m(:,j)',param.covar(:,:,j))*param.w(j);
                    end
                end
                
            case {'stud','studB'}
                %val_pdf=mvtpdf(data-(param.m(:,j)*ones(1,size(data,1)))',param.covar(:,:,j),param.v(j))*param.w(j);%2000=
                if(param.bkg(j))
                    val_pdf     = param.w(j)/prod(sqrt(eig(param_best.covar(:,:,j))));
                else
                    val_pdf     = customStud(data',param.m(:,j), param.covar(:,:,j), param.v(j))*param.w(j);
                end
        end
    catch ME
        disp(ME.message);
    end
end

