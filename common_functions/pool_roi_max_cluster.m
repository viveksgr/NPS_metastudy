function roi_stats = pool_roi_max_cluster(dat, atlas, t_thresh, k_min)
% dat       : fmri_data object (dat.dat is [nVoxels × nSubjects] of t-scores)
% atlas     : canlab atlas (with atlas.labels and atlas.volInfo)
% t_thresh  : voxel‐wise t cutoff
% k_min     : minimum suprathreshold voxels in an ROI
%
% roi_stats : [nSubjects × nROIs] matrix of ROI‐level max-t (NaN if <k_min)

nVox   = size(dat.dat,1);
nSubj  = size(dat.dat,2);
nROIs  = numel(atlas.labels);
roi_stats = nan(nSubj, nROIs);

% Pre‐compute, for each ROI, the linear indices of voxels in the dat.dat mask
for r = 1:nROIs
    % grab the ROI mask as a vector of 0/1 over *the same* in-mask voxels
    roi_img = get_wh_image(atlas, r);       % image_vector: atlas ROI r
    roi_idx = find( roi_img.dat );          % indices into dat.dat’s rows
    % roi_vox_idx{r} = roi_idx;               % store for reuse
end

% Loop subjects
for i = 1:nSubj
    subj_vec = dat.dat(:,i);  % t-value at each voxel for subject i
    
    for r = 1:nROIs
        vals = subj_vec( roi_vox_idx{r} );    % t’s in ROI r
        supr = vals(vals>t_thresh);          % those above threshold
        
        if numel(supr) >= k_min
            roi_stats(i,r) = max(supr);
        else
            roi_stats(i,r) = NaN;             % or zero, your choice
        end
    end
end
end
