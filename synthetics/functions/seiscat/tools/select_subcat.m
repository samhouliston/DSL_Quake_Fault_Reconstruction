function cat = select_subcat(cat, idx)
% Extract sub cat with only selected entries. For fields with neq entries,
% extract all entries listed in <idx>. All other fields leave unchanged
neq  = numel(cat.n);

% If second input argument is logical vector instead of list of indices,
% compute List of indices
if islogical(idx) && size(idx,1)==neq
    idx = find(idx);
end

fdnames  = fieldnames(cat);
nentries = structfun(@(z) size(z,1),cat);
fdnames  = fdnames(nentries==neq);
nf       = numel(fdnames);
for ifd = 1:nf
    
    thisField       = fdnames{ifd};
    cat.(thisField) = cat.(thisField)(idx,:);
end

% Delete ref fields if they exist, since they are no longer correct if the 
% catalogue order is changed
try
    cat = rmfield(cat,'ref');
    cat = rmfield(cat,'ref_rh');
    cat = rmfield(cat,'ref_re');
    cat = rmfield(cat,'ref_azi');

catch
    1;
end