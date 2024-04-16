function clustIDX = mixture2clustIDX(all_points_mat,param)
NUM_P    = size(all_points_mat,1); 
nKernels = numel(param.w);
all_points_matT     = all_points_mat';
prob_mat = zeros(NUM_P,nKernels);
parfor j=1:nKernels
    prob_mat(:,j)    = mvncdf_NOmonte(all_points_mat,param,j,all_points_matT);
end
[~,clustIDX]    = max(prob_mat,[],2);