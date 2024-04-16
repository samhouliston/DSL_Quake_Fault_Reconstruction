function [param_best, best_iCond, best_out_Props]=hierClustSpheres_Fast(all_points_mat,ker,DO_PLOT)
NUM_P           = size(all_points_mat,1);

CUT_AT_DROP     = 1; % Stop when a the first drop off of number of final clusters
MIN_PTS_PLANE   = 4; % Minimum numbers per cluster, smaller ones are killed
MAX_K           = NUM_P/3;%(MIN_PTS_PLANE-1);      %maximum number of clusters
MIN_K           = MAX_K*0.3;                  %minium number of clusters
STEP            = 1;                        %step
SCAN_DIR        = 1; %forward min --> max
%SCAN_DIR        = -1; %backward max --> min
if(nargin<3)
    DO_PLOT     = 0;
end
    if(nargin==2)
        kernel_vec  =ker;
        num_T       = 1;
        CUT_AT_DROP     = 0;
        DO_BBOX     = 1;
    else
        num_T   = round((MAX_K-MIN_K)/STEP);   
        kernel_vec = round(linspace(MIN_K,MAX_K,num_T));
        kernel_vec = unique(kernel_vec);
        num_T      = numel(kernel_vec);
        DO_BBOX     = 0;
    end
    
    if(SCAN_DIR==-1)
        kernel_vec          = kernel_vec(end:-1:1);
    end
    final_kernel_vec   = zeros(1,numel(kernel_vec));

    
    Z = linkage(all_points_mat,'ward','euclidean');
    
    
    for p=1:numel(kernel_vec)
        
        clust_n     = kernel_vec(p);
        clust_id    = cluster(Z,'cutoff',Z(end-(clust_n-1),3),'criterion','distance');
        %clust_id    = cluster(Z,'maxclust',clust_n);
        histVec     = hist(clust_id,1:max(clust_id));
        validVec    = histVec>=MIN_PTS_PLANE;
        NUM_C       = sum(validVec);
        clustCID    = find(validVec);
        bkgRatio    = sum(histVec(~validVec))/NUM_P;
        bkgON       = bkgRatio>0;
        
        param.m         = zeros(3,NUM_C+bkgON);
        param.w         = zeros(1,NUM_C+bkgON);
        param.bkg       = false(1,NUM_C+bkgON); %background noise or not?
        param.sphR      = zeros(1,NUM_C+bkgON);
        param.covar     = zeros(3,3,NUM_C+bkgON);
        param.bbox      = zeros(8,3,NUM_C+bkgON); %bounding box vertices
        param.sig       = zeros(3,NUM_C+bkgON);
        param.mix_type  = 'gaus';
        param.krn_type  = ones(1,NUM_C+bkgON)*2; %1 sphere, 2 gaus, 3 stud;

        for i=1:NUM_C
            tempIND     = clust_id==clustCID(i);    
            tempNUM     = sum(tempIND);
            
            c_points_mat    = all_points_mat(tempIND,:);
            c_mean          = mean(c_points_mat,1);
            sq_dist         = sum(bsxfun(@minus,c_points_mat,c_mean).^2,2);
            c_var           = mean(sq_dist);
            param.m(:,i)    = c_mean';
            param.w(i)      = tempNUM/NUM_P;
            if(tempNUM>=4)
                param.covar(:,:,i)= cov(c_points_mat);
            else
                param.covar(:,:,i)= eye(3)*c_var;       
            end
            param.sig(:,i)       = sqrt(diag(param.covar(:,:,i)));
            
            %For a ball with uniform PDF
            %{
            sphR            = 1.291 * sqrt(c_var); %var=3/5 R^2
            uni_pdf         = 1/((4/3)*pi*(sphR)^3);
            uni_sph_pdf(tempIND,i) = uni_pdf*param.w(i);
            param.covar(:,:,i)  = eye(3)*sphR;
            param.sphR(i)       = sphR;
            %}
        end
        if(bkgON)
            [corPnts, boxCen, covMat]       = minboundboxYK(all_points_mat(:,1),all_points_mat(:,2),all_points_mat(:,3),'v',1);
            param.bkg(NUM_C+1)         = true; 
            param.w(NUM_C+1)           = bkgRatio;
            param.m(:,NUM_C+1)         = boxCen;
            param.covar(:,:,NUM_C+1)   = covMat;
            param.bbox(:,:,NUM_C+1)    = corPnts;
            param.sphR(NUM_C+1)        = 0;
            param.krn_type(NUM_C+1)    = 0;
        end
        %{
        BIC_prev    = nan;
        BIC_ratio   = 1;
        while (BIC_ratio>1e-5 || isnan(BIC_ratio))
            param       = sinpleStepEM_HC(all_points_mat,param,MIN_PTS_PLANE,1,1);
            BIC_ratio   = (BIC_prev-param.BIC)/BIC_prev;
            BIC_prev    = param.BIC;
            disp(['BIC: ' num2str(BIC_prev)]);
        end
        %}
        
        [param_cell,out_Props]  = sinpleStepEM_HC(all_points_mat,param,MIN_PTS_PLANE,DO_BBOX);
        %% Calculate discarted data after EM
        iCond                   = param2iCond(param_cell);
        %{
        [out_Props, param_cell] = eval_mix_model_Fast('gaus',numel(iCond.weightCell),...
                                    all_points_mat',[],...
                                    iCond,1,...
                                    1,0);
        %}
        disp(['Capacity : ' num2str(p) '/' num2str(num_T)...
            '  ' num2str(NUM_C) '/ ' num2str(clust_n)]);
        
        if(out_Props.NumFaults>=max(final_kernel_vec))
            param_best      = param_cell;
            hard_clust      = out_Props.hard_clust;
            best_fin_NUM_C  = out_Props.NumFaults;
            best_NUM_C      = NUM_C;
            best_NLL_Train  = 0;%out_Props.NLL_Train;
            best_iCond      = iCond;
            out_Props.NUMint =NUM_C;
            best_out_Props  = out_Props;
        end                        
        final_kernel_vec(p)     = out_Props.NumFaults;%fin_NUM_C;
        
        if(p>2 && CUT_AT_DROP)
            if(final_kernel_vec(p-1)>final_kernel_vec(p) && ...
                final_kernel_vec(p-2)>final_kernel_vec(p))
                disp(['Output : ' ...
                    '  ' num2str(best_fin_NUM_C) '/ ' num2str(best_NUM_C)]);
                break
                
            end
        end
    end
    
    if(DO_PLOT)
        plot_mix_model(param_best,all_points_mat',[],0);
        title(['NLL' num2str(best_NLL_Train) ' Initial: ' num2str(best_NUM_C) ' Final: ' num2str(best_fin_NUM_C)])
    end
    
    %%
    %{
    if(DO_PLOT)
        figure;
        subplot(2,1,1);
        scatter(kernel_vec,final_kernel_vec);
        initID=find(final_kernel_vec==max(final_kernel_vec),1);
        title(['MAX: '  num2str(max(final_kernel_vec)) ' @ ' num2str(kernel_vec(initID))]);
        grid on;
        xlabel('Initial Clusters');
        ylabel('Final Clusters');
    end
    %}
end