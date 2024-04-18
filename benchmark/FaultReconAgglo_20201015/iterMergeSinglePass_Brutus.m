function param_best=iterMergeSinglePass_Brutus(param_best,Pts_Train,MergeUNR,OmitNMDL,hard_clust,GAIN_MODE)
%First pass: consider only events within candidate merging couples
%Second pass: consider all events within the catalog 

% Y.Kamer 20201015

DO_PLOT_GAINMAT = 0; %plot gain matrix

DO_ITER_PLOT    = 0; %plot iteration results
SAVE_PLOT       = 0; %save iteration results as seprate images

SECOND_PASS     = 1;
if(nargin<6)
    GAIN_MODE   = 'local';
end
if(strcmp(GAIN_MODE,'localOnly'))
    GAIN_MODE       = 'local';
    SECOND_PASS     = 0;
    MINpdf          = eps;%eps(0);
else
    MINpdf          = eps;%eps(0);
end

CHECK_SIGMA  = 1; % Check if bounding boxes are touching, if not dont consider merger

DO_PLOT     = 0;
VERBOSE     = 1;
tic;

FIRST_LOOP = 1; 
NO_PTS      = size(Pts_Train,2);
nKernels    = numel(param_best.w);
gainMat     = zeros(nKernels);

if(VERBOSE)
    disp([datestr(now) '  Calculate initial likelihood']);
end

[init_PDF,init_MDL] = MDLGaus(Pts_Train',param_best,[],MINpdf);

%disp(init_PDF(1:20))
gain_best   = 0;
iterCNT     = 0;
param_best.sig = nan(size(param_best.m));
for i=1:nKernels
    param_best.sig(:,i)= sqrt(diag(param_best.covar(:,:,i))); 
end


while(1)
    iterCNT     = iterCNT+1;
    BRUTUS_DONE = 0;
    %
    is_na = isnan(gainMat);
    disp(['Number of nans in gainMat ' num2str(sum(is_na(:))/2) '/' num2str(size(gainMat,1)*(size(gainMat,1)-1)/2)]);
    [iV,jV]     = ind2sub(size(gainMat),find(triu(ones(size(gainMat)),1)~=0 & ~isnan(gainMat)));
    pair_vec    = [iV jV];
    NUM_pr      = size(pair_vec,1);

    if(VERBOSE)
        disp([datestr(now) '  Clusters >> ' num2str(nKernels) ', MDL: ' num2str(init_MDL) '           ']);
    end
    % renove pairs without overlapping bbox
    if(CHECK_SIGMA)
        TF_vec      = false(NUM_pr,1);
        parfor i=1:size(pair_vec,1)
            if(mod(i,100)==0);progScreen(i);end
            TF_vec(i)=touchBBOX(param_best,iV(i),jV(i));
        end
        
        disp([num2str(sum(TF_vec)) '/'  num2str(NUM_pr) ' pairs have touching bbox']);
    end
    if(isempty(pair_vec))
        pair_vec=ones(0,2);
    end
    if(CHECK_SIGMA)
        idxN = sub2ind(size(gainMat), pair_vec(~TF_vec,1),pair_vec(~TF_vec,2));
        gainMat(idxN)=NaN;
        disp(['Number of nans from bbox check: ' num2str(length(idxN))]);
        idxN = sub2ind(size(gainMat), pair_vec(~TF_vec,2),pair_vec(~TF_vec,1));
        gainMat(idxN)=NaN;
        pair_vec    = pair_vec(TF_vec,:);
    end
    is_na = isnan(gainMat);
    disp(['Number of nans after bbox check: ' num2str(sum(is_na(:))/2)]);
    
    res_vec  = zeros(size(pair_vec,1),1);
    if(VERBOSE)
        disp([num2str(iterCNT) '>> Pair merges to calc: ' num2str(numel(res_vec)) '           ']);
    end
    
    switch GAIN_MODE
        case 'local' %consider only events in the merging clusters
            parfor n=1:size(pair_vec,1)%par
                if(mod(n,100)==0);progScreen(n);end
                c1          = pair_vec(n,1);
                c2          = pair_vec(n,2);
                e1          = hard_clust{c1}.num_e;
                e2          = hard_clust{c2}.num_e;
                dataID      = hard_clust{c1}.ID|hard_clust{c2}.ID;
                prm_m       = mergeGausOne(param_best,c1,c2);
                [~,mrg_MDL] = MDLGaus(Pts_Train(:,dataID)',prm_m,[],MINpdf);
                sep_MDL     = MDLGausSel(Pts_Train(:,dataID)',param_best,[c1 c2],[e1 e2],MINpdf);
                res_vec(n)  = sep_MDL-mrg_MDL;
            end
                        
        case 'global'%consider all events
            %if pair number low work localy, else call Brutus
            iV  = pair_vec(:,1);
            jV  = pair_vec(:,2);
            parfor n=1:size(pair_vec,1)%par
                if(mod(n,100)==0);progScreen(n);end
                c1          = iV(n);
                c2          = jV(n);
                sep_PDF     = MDLGaus(Pts_Train',param_best,[c1 c2]);
                prm_m       = mergeGausOne(param_best,c1,c2,'KEEP_WGT');
                mrg_PDF     = MDLGaus(Pts_Train',prm_m,1);
                sum_PDF     = init_PDF - sep_PDF + mrg_PDF;
                sum_PDF(sum_PDF<=MINpdf) = MINpdf;
                perMDL  = sum(-log(sum_PDF));
                MT      = 1; % 1 for Gaus, 2 for Stud
                NO_KRN  = nKernels-1;
                prm_MDL  = perMDL+0.5*(MT*NO_KRN-1+9*NO_KRN)*log(NO_PTS);
                res_vec(n)  = init_MDL-prm_MDL;
            end
    end
    
    idx = sub2ind(size(gainMat), pair_vec(:,1),pair_vec(:,2));
    gainMat(idx)=res_vec(1:end);
    idx = sub2ind(size(gainMat), pair_vec(:,2),pair_vec(:,1));
    gainMat(idx)=res_vec(1:end);
    if(DO_PLOT_GAINMAT)
        figure;
        gainMat(isnan(gainMat))=0;
        imagesc(tril(gainMat));
        daspect([1 1 1]);
        cmap     = cbrewer('div','RdYlGn',11);
        cmap(6,:)=[1 1 1];
        colormap(cmap);
        colorbar;
        caxis([-20 20])
    end
    
    gain_best   = max(gainMat(1:end));
    if(~(gain_best>0))
        switch GAIN_MODE
            case 'local' % finished local optimization now do global
                if(SECOND_PASS)
                    FIRST_LOOP  = 1;
                    GAIN_MODE   = 'global';
                    disp('>> Switching to Global MLD')
                    continue;
                else
                    break;
                end
            case 'global'% global is finished nothing more to do
                break
        end
    end
    [mi,mj]     = find(gainMat==gain_best,1);
    vi=mi;
    vj=mj;
    inC = 1;
    ID_posMDL   = res_vec>0;
    
    disp([num2str(sum(ID_posMDL)) '/' num2str(length(ID_posMDL)) ' pairs with positive gain'])
    pair_posMDL = pair_vec(ID_posMDL,:);
    mdl_posMDL  = res_vec(ID_posMDL);
    if(OmitNMDL)
        pair_negMDL = pair_vec(~ID_posMDL,:);
        idxN = sub2ind(size(gainMat), pair_negMDL(:,1),pair_negMDL(:,2));
        disp(['Nans from negative pairs: ' num2str(length(idxN))])
        gainMat(idxN)=NaN;
        idxN = sub2ind(size(gainMat), pair_negMDL(:,2),pair_negMDL(:,1));
        gainMat(idxN)=NaN;
    end
    is_na = isnan(gainMat);
    disp(['Number of nans after negative check: ' num2str(sum(is_na(:))/2)]);
    %% Find & exclude cluster which could be merged with the best pair
    while(inC<MergeUNR+1)%
        ex_clust    = find(max(gainMat(mi,:),gainMat(mj,:))>0);
        ID_del      = ismember(pair_posMDL,ex_clust);
        ID_del      = max(ID_del,[],2);
        pair_posMDL = pair_posMDL(~ID_del,:); 
        mdl_posMDL  = mdl_posMDL(~ID_del);
        if(isempty(pair_posMDL))
            break
        end
        tmp     = pair_posMDL(find(mdl_posMDL==max(mdl_posMDL),1),:);
        mi      = tmp(1);
        mj      = tmp(2);
        vi      =[vi mi];
        vj      =[vj mj];
        inC     =inC+1;
    end

    sep_PDF_T=0;
    mrg_PDF_T=0;
    N_MRG    = numel(vi);
    %{
    if(any(param_best.bkg([vi vj])))
        ('Background merging!')
    end
    %}
    disp([num2str(N_MRG) ' pairs merged'])
    for i=1:N_MRG
        sep_PDF     = MDLGaus(Pts_Train',param_best,[vi(i) vj(i)]);
        prm_m       = mergeGausOne(param_best,vi(i),vj(i),'KEEP_WGT');
        mrg_PDF     = MDLGaus(Pts_Train',prm_m,1);
        sep_PDF_T   = sep_PDF_T+sep_PDF;
        mrg_PDF_T   = mrg_PDF_T+mrg_PDF;
    end
    %Check if MDL is going down, if not=> break
    T_init_PDF      = init_PDF - sep_PDF_T + mrg_PDF_T;
    T_init_PDF(T_init_PDF<=MINpdf)    = MINpdf;
    T_init_MDL      = sum(-log(T_init_PDF));
    MT              = 1; % 1 for Gaus, 2 for Stud
    NO_KRN          = nKernels-N_MRG;
    T_init_MDL      = T_init_MDL+0.5*(MT*NO_KRN-1+9*NO_KRN)*log(NO_PTS);
    if(T_init_MDL==inf)
        break;
        continue
    else
        init_MDL    = T_init_MDL;
        init_PDF    = T_init_PDF;
        param_best  = mergeGaus(param_best,vi,vj);
        if(strcmp(GAIN_MODE,'local'))
            hard_clust  = mergeHard(hard_clust,vi,vj);
        end
    end
        

    
    if(FIRST_LOOP);FIRST_LOOP=0;end
    
    is_na = isnan(gainMat);
    disp(['Nans right before removing merged entries ' num2str(sum(is_na(:))/2)]);
    % Remove the merged clusters from the matrix
    %B = gainMat;
    gainMat([vi vj],:)=[];
    gainMat(:,[vi vj])=[];
    % Add a row-column for the new cluster
    [r,~]       = size(gainMat);
    nKernels    = r+inC;
    
    is_na = isnan(gainMat);
    disp(['Nans right after removing merged entries ' num2str(sum(is_na(:))/2)]);

    if(isempty(gainMat))
        gainMat =   zeros(nKernels);
    else
        gainMat = [gainMat zeros(r, inC)];
        gainMat = [gainMat; zeros(inC, nKernels)];
        %gainMat(r:nKernels,:)=0;
        %gainMat(:,r:nKernels)=0;
    end
    
    is_na = isnan(gainMat);
    disp(['Nans after adding new entries for new kernels ' num2str(sum(is_na(:))/2)]);
    
    %% Plot the merged network
    if(DO_ITER_PLOT)
        plot_mix_model(param_best,[],[],[],0,0);
        title([num2str(nKernels) ' Clusters']);
        pause(1);
        if(SAVE_PLOT)
            axis tight;
            %Landers Plots
            %{
            ylim([0 100]);
            xlim([5 55]);
            zlim([0 17]);
            set(gca,'ZDir','reverse');
            set(gca,'YDir','normal');
            view(65,36)
            set(gca,'xtick',[],'ytick',[],'ztick',[])
            grid off;
            box on;
            %}
            %Synth 5 faults Plots
            %
            axis([-10 58 -20 25 -5 12])
            grid on;
            box on;
            view(25,50);
            %}
            set(gcf,'InvertHardcopy','off','color',[1 1 1]*0.8)
            print('-dpng','-r250',['Merge_Iter_' num2str(iterCNT,'%03d') '.png']);
            close;
        end
    end
end
eTm=toc;

param_best.mdl  = init_MDL;
nKernels = numel(param_best.w);
disp(['Merged atoms into: ' num2str(nKernels) ', time= ' num2str(eTm)]);

    if(DO_PLOT)
        plot_mix_model(param_best,[],[],0);
        title([num2str(nKernels) ' Cluster , MDL = ' num2str(init_MDL)]);
    end
end