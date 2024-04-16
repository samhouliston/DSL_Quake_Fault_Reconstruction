function perMDL = MDLGausSel(all_points_mat,param,vecIND,numE,MINpdf)
    NO_PTS      = size(all_points_mat,1);
    nKernels    = numel(vecIND);
    per_pdf     = zeros(NO_PTS,1);
    param.w(vecIND)=numE/sum(numE);
    all_points_matT = all_points_mat';
    for j=1:nKernels
        tmp_pro       = mvncdf_NOmonte(all_points_mat,param,vecIND(j),all_points_matT);
        per_pdf(:)    = per_pdf(:) + tmp_pro;
    end
    %Apply zero correction only if MDL is to be calculated
    zrID    = per_pdf<MINpdf;
    if(sum(zrID))
        disp('###################<MDLGausSel>  MINpdf Detected!')
        per_pdf(zrID)   = MINpdf;
        snd = load('gong.mat');
        sound(snd.y,10*snd.Fs);
    end
    perMDL  = sum(-log(per_pdf));
    MT      = 1; % 1 for Gaus, 2 for Stud 
    perMDL  = perMDL+0.5*(MT*nKernels-1+9*nKernels)*log(NO_PTS);
end
