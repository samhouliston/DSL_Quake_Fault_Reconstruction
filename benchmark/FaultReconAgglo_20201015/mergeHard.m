function hard_clust=mergeHard(hard_clust,vi,vj)
    for k=1:numel(vi)
        i=vi(k);
        j=vj(k);
        ID      = hard_clust{i}.ID|hard_clust{j}.ID;
        num_e   = hard_clust{i}.num_e + hard_clust{j}.num_e;
        
        numK   = numel(hard_clust);
        hard_clust{numK+1}.ID       = ID;
        hard_clust{numK+1}.num_e    = num_e;
    end
    hard_clust([vi vj])           = [];
end