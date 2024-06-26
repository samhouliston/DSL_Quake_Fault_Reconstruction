% This script calls the agglomorative fault reconstruction method on a
% given dataset. The method is described in the following paper:
% Kamer Y., Ouillon G., Sornette D. (2020) "Fault Network Reconstruction using Agglomerative Clustering:
% Applications to South Californian Seismicity"

% Y.Kamer 20201015

% Changes made:
% - changed MINpdf to eps from eps(0)
% - buffered and rounded values in minboundboxYK

GAIN_MODE       ='global';  % Merger gain mode, set to 'global' or 'localOnly'
CHUNK_NO        = 4;        % Split the data into multiple chunks, set to 1 for no split 
DELETE_BKG      = 0;        % [!Experimental!] If set to 1: Deletes background clusters and points associated with them, set to 0 to keep the original data 

PLOT_CHUNK      = 0;        % Plot results of individual chunks 
YDIR            = 'normal'; % Direction of y axis, set to 'normal' or 'reverse'

%RUN_CASE = 'SimpleGauss';
RUN_CASE = 'Landers';
%RUN_CASE = 'Synth_5Faults';

%% Add your run case to the following block
switch RUN_CASE
    case 'Landers'
        % load Landers hypocenters
        ld              = load('landers.mat');
        all_points_mat  = ld.xyz_mat;
    case 'Synth_5Faults'
        % Generate synthetic faults
        SPC_FAC         = 2;    % Spacing factor, higher values -> faults more seperated
        DENS_FAC        = 0.5;  % Density factor, higher values -> more samples
        MODE_ELIP       = 0;    % 0: Rectangular planes, 1: Gaussian ellipses 
        NUM_FAULTS      = 5;    % Number of faults
        NOISE_PRC       = 20;   % Background noise percentage
        RND_SEED        = 5;    % Random number generator seed
        DO_PLOT         = 1;    % 0: no plot, 1: colored plot, 2: unmarked plot

        catNUM          = numel(NUM_FAULTS);
        cellPTS         = cell(catNUM,1);
        figure;
        title('');
        all_points_mat  = syth_fault_events(NUM_FAULTS, RND_SEED, NOISE_PRC, DENS_FAC,...
                            SPC_FAC, MODE_ELIP,DO_PLOT);
        view(25,50);
    case 'SimpleGauss'
        ld = load('ex_dat.mat');
        all_points_mat = ld.xyz;
end

%%
if(PLOT_CHUNK)
    figure;
    plot3(all_points_mat(:,1),all_points_mat(:,2),all_points_mat(:,3),'.r');
    daspect([1 1 1])
    set(gca,'ydir',YDIR)
    set(gca,'zdir','reverse')
    axis tight
    view(2);
    set(gca,'color',[1 1 1]*0.2)
    set(gca,'userdata','Main')
    h(1)    = gca;
    hMain   = findobj('userdata','Main');
end


%% Build Hierarchical tree; cut into a few big clusters
NUM_PTS     = size(all_points_mat,1);
param_AllBKG= [];
Z           = linkage(all_points_mat,'ward','euclidean');
clust_id    = cluster(Z,'cutoff',Z(end-(CHUNK_NO-1),3),'criterion','distance');
hist_vec    = hist(clust_id,1:max(clust_id));
chunk_vec   = find(hist_vec);
hist_vec    = hist_vec(chunk_vec);
figure;
bar(sort(hist_vec))
title('Number of points per cluster');
%%
num_chunk   = numel(chunk_vec);
cmap        = hsv(5);
cmap        = cmap(randperm(num_chunk),:);
[~,IX]      = sort(hist_vec);
figure;
for i=1:num_chunk
    tempIND         = clust_id==chunk_vec(IX(i));
    chunk_mat       = all_points_mat(tempIND,:);
    plot(chunk_mat(:,1),chunk_mat(:,2),'.k',...
                'markersize',3,...
                'Color',cmap(i,:));
    hold on;
end
axis tight;
daspect([1 1 1]);
set(gca,'color',[1 1 1]*0.3);
set(gcf,'InvertHardcopy','off','color',[1 1 1]*1)
title ([num2str(num_chunk) ' Subsets'],'fontsize',14)
%%
if(PLOT_CHUNK)
    fig_c   = figure;
end
indChunk    = 0;
for i=1:num_chunk
    
    tempIND         = clust_id==chunk_vec(i);    
    tempNUM         = sum(tempIND);
    chunk_mat       = all_points_mat(tempIND,:);
    if(tempNUM<5)
        continue;
    end
    if(PLOT_CHUNK)
        figure(fig_c)
        hold on;
        plot3(chunk_mat(:,1),chunk_mat(:,2),chunk_mat(:,3),'.k',...
                'Color',cmap(i,:));
        
        hA = findall(gcf,'type','axes');
        set(hA,'Xlim',get(hMain(1),'Xlim'),'Ylim',get(hMain(1),'Ylim'),'Zlim',get(hMain(1),'Zlim'));
        daspect([1 1 1]);
        title([ 'Chunk: ' num2str(i) '/' num2str(num_chunk)]);
        set(gca,'color',[1 1 1]*0.2)
        
        set(gca,'ydir',YDIR)
        set(gca,'zdir','reverse')
        axis tight
    end
    %% Atomize each chunk
    numMAXC              = hierClustSpheres_JustNum(chunk_mat);
    [param_best,iCond, out_Props]   = hierClustSpheres_Fast(chunk_mat,numMAXC);
    
    hard_clust      = out_Props.hard_clust;
    idBKG           = find(param_best.bkg);

    if(DELETE_BKG && ~isempty(idBKG))
        if(isempty(param_AllBKG))
            param_AllBKG.m          = param_best.m(:,idBKG);
            param_AllBKG.w          = param_best.w(idBKG)*(tempNUM/NUM_PTS);
            param_AllBKG.bkg        = param_best.bkg(idBKG);
            param_AllBKG.sphR       = param_best.sphR(idBKG);
            param_AllBKG.covar      = param_best.covar(:,:,idBKG);
            param_AllBKG.bbox       = param_best.bbox(:,:,idBKG);
            param_AllBKG.sig        = param_best.sig(:,idBKG);
            param_AllBKG.krn_type   = param_best.krn_type(idBKG);
        else
            curNUM  = numel(param_AllBKG.w)+1;
            
            param_AllBKG.m(:,curNUM)        = param_best.m(:,idBKG);
            param_AllBKG.w(curNUM)          = param_best.w(idBKG)*(tempNUM/NUM_PTS);
            param_AllBKG.bkg(curNUM)        = param_best.bkg(idBKG);
            param_AllBKG.sphR(curNUM)       = param_best.sphR(idBKG);
            param_AllBKG.covar(:,:,curNUM)  = param_best.covar(:,:,idBKG);
            param_AllBKG.bbox(:,:,curNUM)   = param_best.bbox(:,:,idBKG);
            param_AllBKG.sig(:,curNUM)      = param_best.sig(:,idBKG);
            param_AllBKG.krn_type(curNUM)   = param_best.krn_type(idBKG);
        end
        param_best.m(:,idBKG)            = [];
        param_best.w(idBKG)              = [];
        param_best.bkg(idBKG)            = [];
        param_best.sphR(idBKG)           = [];
        param_best.covar(:,:,idBKG)      = [];
        param_best.bbox(:,:,idBKG)       = [];
        param_best.sig(:,idBKG)          = [];
        param_best.krn_type(idBKG)       = [];
        delPTS  = false(size(hard_clust{1}.ID));
        for d=1:numel(idBKG)
            delPTS = delPTS | hard_clust{idBKG(d)}.ID;
        end
        hard_clust(idBKG)=[];
        chunk_mat(delPTS,:)=[];
        for d=1:numel(hard_clust)
            hard_clust{d}.ID(delPTS)=[];
        end
    end
    
    save('params_preMerge.mat','-struct','param_best')

    param_Inc       = iterMergeSinglePass_Brutus(param_best,chunk_mat',inf,true,hard_clust,GAIN_MODE);
    disp(['Chunk ID: '  num2str(i) ' ^^']);

    if(PLOT_CHUNK)
        %
        set(gca,'color',[1 1 1]*0.2)
        hA=findall(gcf,'type','axes');
        set(hA,'Xlim',get(hMain(1),'Xlim'),'Ylim',get(hMain(1),'Ylim'),'Zlim',get(hMain(1),'Zlim'))
        set(gca,'ydir',YDIR)
        set(gca,'zdir','reverse')
        axis tight
        %}
    end
    
    [param_Inc,out_Props]   = sinpleStepEM_HC(chunk_mat,param_Inc,1,1);
    hard_clust              = out_Props.hard_clust;
    % Adjust the weights of each chunk
    param_Inc.w = param_Inc.w*(size(chunk_mat,1)/NUM_PTS);
    
    if(i==1)
        param_All       = param_Inc;
        hard_clust_All  = hard_clust;
        chunk_mat_All   = chunk_mat;
    else
        param_All       = paramJoin(param_All,param_Inc);
        hard_clust_All  = hardJoin(hard_clust_All,hard_clust);
        chunk_mat_All   = [chunk_mat_All;chunk_mat];
    end
    
    %disp(param_All.m)
    %disp(param_All.covar)
    %disp(param_All.w)
    %disp(param_All.bbox)
end

if(PLOT_CHUNK)
    figure(fig_c)
    title([num2str(num_chunk) ' Chunks']);
end
        
%% Merge all the accumulated chunk
if(num_chunk~=0)
    
    %% Plot initial input
    %{
    plot_mix_model(param_All,[],[],0,1);
    title(['INITIAL: ' num2str(numel(param_All.w)) ' clusters']);
    view(2);
    whitebg([1 1 1]*1)
    %}
    %%
    switch GAIN_MODE
        case 'localOnly'
            param_All_FIN_LOC           = iterMergeSinglePass_Brutus(param_All,chunk_mat_All',inf,true,hard_clust_All,'localOnly');
            % Plot LOCAL
            plot_mix_model(param_All_FIN_LOC,[],[],[],1,0);
            title(['LOCAL: ' num2str(numel(param_All_FIN_LOC.w)) ' clusters, BKG = ' num2str(sum(param_All_FIN_LOC.w(param_All_FIN_LOC.bkg)),'%.2f')],'fontsize',14);
        case 'global'
            param_All_FIN_GLO           = iterMergeSinglePass_Brutus(param_All,chunk_mat_All',inf,true,[],'global');
            save('params_postMerge.mat','-struct','param_All_FIN_GLO');
            % Plot GLOBAL
            disp(param_All_FIN_GLO.w(param_All_FIN_GLO.bkg))
            plot_mix_model(param_All_FIN_GLO,[],[],[],1,0);
            title(['GLOBAL: ' num2str(numel(param_All_FIN_GLO.w)) ' clusters, BKG = ' num2str(sum(param_All_FIN_GLO.w(param_All_FIN_GLO.bkg)),'%.2f')],'fontsize',14);        
    end
    %view(2);
    %whitebg([1 1 1]*1)
    
    
    
end
