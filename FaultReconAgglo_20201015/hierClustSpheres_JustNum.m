function best_NUM_C = hierClustSpheres_JustNum(all_points_mat,ker)
NUM_P           = size(all_points_mat,1);

CUT_AT_DROP     = 0; % Stop when a the first drop off of number of final clusters
MIN_PTS_PLANE   = 4; % Minimum numbers per cluster, smaller ones are killed
MAX_K           = NUM_P/3;%(MIN_PTS_PLANE-1);      %maximum number of clusters
MIN_K           = MAX_K*0.3;                  %minium number of clusters
STEP            = 1;                        %step
SCAN_DIR        = 1; %forward min --> max
%SCAN_DIR        = -1; %backward max --> min
DO_PLOT     = 0;
    if(nargin==2)
        kernel_vec  =ker;
        num_T       = 1;
        CUT_AT_DROP     = 0;
        DO_BBOX     = 1;
    else
        num_T   = round((MAX_K-MIN_K)/STEP);   
        disp(num_T);
        kernel_vec = round(linspace(MIN_K,MAX_K,num_T));
        kernel_vec = unique(kernel_vec);
        num_T      = numel(kernel_vec);
        DO_BBOX     = 0;
    end
    
    if(SCAN_DIR==-1)
        kernel_vec          = kernel_vec(end:-1:1);
    end
    final_kernel_vec    = zeros(1,numel(kernel_vec));
    final_bkg           = zeros(1,numel(kernel_vec));
    
    Z = linkage(all_points_mat,'ward','euclidean');
    
    
    for p=1:numel(kernel_vec)
        
        clust_n     =kernel_vec(p);
        clust_id    = cluster(Z,'cutoff',Z(end-(clust_n-1),3),'criterion','distance');
        %clust_id    = cluster(Z,'maxclust',clust_n);
        NUM_C       = max(clust_id);
        histVEC     = hist(clust_id,1:NUM_C);
        valID       = histVEC>=MIN_PTS_PLANE;
        fin_NUM_C   = sum(valID);
        final_bkg(p)= sum(histVEC(~valID))/NUM_P;
        %{
        disp(['Capacity : ' num2str(p) '/' num2str(num_T)...
            '  ' num2str(fin_NUM_C) '/ ' num2str(NUM_C)]);
        %}
        if(fin_NUM_C>=max(final_kernel_vec))
            best_fin_NUM_C  = fin_NUM_C;
            best_NUM_C      = NUM_C;
        end                        
        final_kernel_vec(p)     = fin_NUM_C;
        
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
        figure;
        %subplot(2,1,1)
        ttl = [ 'Capacity = ' num2str(max(final_kernel_vec)) ' @' num2str(best_NUM_C) ];
        plot(kernel_vec,final_kernel_vec,'.k','displayname',ttl);
        xlabel('Initial Clusters');
        ylabel('Final Clusters');
        
        %{
        subplot(2,1,2)
        plot(kernel_vec,final_bkg,'.k');
        xlabel('Initial Clusters');
        ylabel('Background noise ratio');
        %}
    end
end