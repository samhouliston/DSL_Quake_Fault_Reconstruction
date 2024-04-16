function [param, out_Props,chrome] = sinpleStepEM_HC(all_points_mat,param,MIN_PTS_PLANE,DO_BBOX,DO_EXP,IN_chrome)
%single step expectation maximization hard clustering
NUM_P    = size(all_points_mat,1); 
nKernels = numel(param.w);
MINpdf   = eps; %eps(0);
all_points_matT     = all_points_mat';
if(nargin<5)
    DO_EXP=0;
end
if(nargin==6)
    chrome = IN_chrome;
else
    prob_mat = zeros(NUM_P,nKernels);
    parfor j=1:nKernels
        prob_mat(:,j)    = mvncdf_NOmonte(all_points_mat,param,j,all_points_matT);
    end
    %% Assign points to their mixture kernels and plot planes for each cluster
    [~,chrome]    = max(prob_mat,[],2);
   
end

sum_e           = 0;
IDdel           = true(1,nKernels);
NUM_vec         = nan(1,nKernels); %number points per kernel
hard_clust      = cell(1,nKernels);

for i=1:nKernels
    IDX_temp        = chrome==i;
    tmpNumE         = sum(IDX_temp);
    NUM_vec(i)      = tmpNumE;
    if(tmpNumE>=MIN_PTS_PLANE || (tmpNumE>0 && param.bkg(i)))
        IDdel(i)                = false;
        hard_clust{i}.ID        = IDX_temp;
        hard_clust{i}.num_e     = tmpNumE;
        
        if(DO_EXP && ~param.bkg(i))
            newCovar                    = cov(all_points_mat(IDX_temp,:));
            if(isinf(cond(newCovar)))
                disp('Fuu Singuu!')
                continue;
            end
            param.covar(:,:,i)          = newCovar;
            param.m(:,i)                = mean(all_points_mat(IDX_temp,:));
        end
        if(tmpNumE>3 && DO_BBOX && ~param.bkg(i)) %Don't update bbox for background
            %
            bboxDATA            = minboundboxYK(all_points_mat(IDX_temp,1),...
                                                all_points_mat(IDX_temp,2),...
                                                all_points_mat(IDX_temp,3),...
                                                'v',1);
            bboxMODEL           = DrawCuboidYK(param.m(:,i),param.covar(:,:,i)*12)';
            vertJ               = [bboxDATA; bboxMODEL]';
            corPts              = minboundboxYK(vertJ(1,:),vertJ(2,:),vertJ(3,:),'v',1);      
            param.bbox(:,:,i)   = corPts;
        end
        param.w(i)                  = tmpNumE;
        sum_e                       = sum_e + tmpNumE;
    end
end
del_vec                    = find(IDdel);
hard_clust(del_vec)        = [];
NUM_vec(del_vec)           = []; 
param.w(del_vec)           = [];
param.bkg(del_vec)         = [];
param.m(:,del_vec)         = [];
param.covar(:,:,del_vec)   = [];
param.bbox(:,:,del_vec)    = [];
param.krn_type(del_vec)    = [];
param.sphR(del_vec)        = [];
param.w=param.w/NUM_P;

prob_mat = zeros(NUM_P,nKernels);
parfor j=1:nKernels
    prob_mat(:,j)    = mvncdf_NOmonte(all_points_mat,param,j,all_points_matT);
end
[~,chrome]      = max(prob_mat,[],2);
per_pdf         = sum(prob_mat,2);
zrID                = per_pdf<MINpdf;
if(sum(zrID))
    disp('##########################  MINpdf Detected!')
    per_pdf(zrID)   = MINpdf;
    snd = load('gong.mat');
    sound(snd.y,10*snd.Fs);
end    
param.BIC                   = sum(-log(per_pdf))+0.5*(nKernels-1+9*nKernels)*log(NUM_P);

out_Props.NumFaults     = numel(param.w);
out_Props.NumVec        = NUM_vec;
out_Props.NumPoints     = NUM_P;
out_Props.hard_clust    = hard_clust;
