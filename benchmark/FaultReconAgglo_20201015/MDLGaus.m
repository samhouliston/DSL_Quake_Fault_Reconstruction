function [per_pdf, perMDL] = MDLGaus(all_points_mat,param,vecKern,MINpdf)
%all_point_mat:     datapoints for likelihood calculation 
%param:             structure of all kernels
%vecKern:           if provided calc likelihood only for these kernels and set MDL to zero

    NO_PTS              = size(all_points_mat,1);
    if(isempty(vecKern))
        nKernels        = numel(param.w);
        vecKern         = 1:nKernels;
        FLAG_SUBSET     = 0;
    else
        nKernels        = numel(vecKern);
        FLAG_SUBSET     = 1;
        perMDL          = nan;
    end
    
    per_pdf     = zeros(NO_PTS,1);
    %disp([datestr(now) '  Likelihood of ' num2str(nKernels) ' clusters']);
    %disp('   ');
    sep_pdf = zeros(NO_PTS, nKernels);
    if(FLAG_SUBSET)
        for j=vecKern
            %progScreen(j);
            tmp_pro     = mvncdf_NOmonte(all_points_mat,param,j,all_points_mat');
            per_pdf     = per_pdf + tmp_pro;
        end
    else
        parfor j=vecKern % if its not a subset then it is sequential
            %progScreen(j);
            tmp_pro     = mvncdf_NOmonte(all_points_mat,param,j,all_points_mat');
            per_pdf     = per_pdf + tmp_pro;
            sep_pdf(:,j) = tmp_pro;
        end
    end
    
    if(nargout==2)
        %Apply zero correction only if MDL is to be calculated
        zrID    = per_pdf<MINpdf;
        if(sum(zrID))
            disp('#################<MDLGaus>  MINpdf Detected!')
            per_pdf(zrID)   = MINpdf;
            snd = load('gong.mat');
            sound(snd.y,10*snd.Fs);
        end


        perMDL  = sum(-log(per_pdf));
        MT      = 1; % 1 for Gaus, 2 for Stud 
        perMDL  = perMDL+0.5*(MT*nKernels-1+9*nKernels)*log(NO_PTS);
    end
end
