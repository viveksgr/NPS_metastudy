function cl = extract_roi_tstats_max(dat, atlas, t_thresh, k_min)
% EXTRACT_ROI_TSTATS_MAX  ROI‐wise max-t cluster summary + group t-test
%
%   cl = extract_roi_tstats_max(dat, atlas, t_thresh, k_min)
%
% Inputs:
%   - dat       : fmri_data object with dat.dat = [nVoxels×nSubjects] of first-level t-maps
%   - atlas     : canlab atlas object
%   - t_thresh  : voxel-wise t threshold (e.g. 3.1)
%   - k_min     : minimum number of suprathreshold voxels in an ROI
%
% Output:
%   - cl(i).label : ROI name
%   - cl(i).t     : group-level t-statistic for ROI i
%   - cl(i).p     : two-tailed p-value
%   - cl(i).df    : degrees of freedom (nEffectiveSubjects−1)

    % pull ROI masks from atlas
    regions = atlas2region(atlas);
    nR = numel(regions);
    nSubj = size(dat.dat,2);

    % precompute voxel indices for each ROI
    for r=1:nR
        % get a binary image_vector for region r
        iv = get_wh_image(atlas, r);  
        % find which rows of dat.dat correspond to that ROI
        roi_idx{r} = find(iv.dat);      
        cl(r).label = regions(r).shorttitle;
    end

    % for each ROI, compute each subject's ROI‐max score
    roi_scores = nan(nSubj, nR);
    for r=1:nR
      idx = roi_idx{r};
      for i=1:nSubj
        vox_t = dat.dat(idx,i);
        supr = vox_t(vox_t > t_thresh);
        if numel(supr) >= k_min
          roi_scores(i,r) = max(supr);
        else
          roi_scores(i,r) = NaN;
        end
      end
    end

    % now run a one-sample t-test across subjects for each ROI
    for r=1:nR
      y = roi_scores(:,r);
      y = y(~isnan(y));           % drop subjects with no cluster
      cl(r).df    = numel(y)-1;
      if cl(r).df > 0
        [~,p,~,stats] = ttest(y, 0);
        cl(r).t  = stats.tstat;
        cl(r).p  = p;
      else
        cl(r).t = NaN;
        cl(r).p = NaN;
      end
    end
end
