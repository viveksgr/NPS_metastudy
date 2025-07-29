function  parc_data = apply_parc_contrasiobj(fmri_objs,atlas,f_correct)

% --- 0. Inputs -----------------------------------------------------------
%   fmri_objs   : 1×11 cell of fmri_data objects
%   f_correct   : 1×11 cell, each a vector of ROI‐indices lost for that study
%   atlas       : canlab atlas struct (with atlas.labels, numel = K)

nStudies = numel(fmri_objs);
nROIs    = numel(atlas.labels);

% Preallocate
parc_data = cell(1, nStudies);

for s = 1:nStudies
    fprintf('Study: %02d\n',s)
    dat = fmri_objs{s};
    
    % 1) extract whatever surviving ROIs you get back
    % roi_dats = extract_roi_averages(dat, atlas);
    roi_dats = extract_roi_averages_maxt(dat, atlas);
    % roi_dats3 = extract_roi_averages_maxt(dat, atlas);

    nFound = length(roi_dats);
    nSubj = length(roi_dats(1).dat);
    roi_means = zeros( nSubj,nFound);
    for nn=1:nFound; roi_means(:,nn)=roi_dats(nn).dat; end
    
    % 2) figure out which ROIs should be missing here
    missing = f_correct{s};              % e.g. [451 455 …] or []
    valid   = setdiff(1:nROIs, missing); % the ROI‐slots you *do* have data for
    
    % sanity check
    if numel(valid) ~= nFound
        error('Study %d: expected %d surviving ROIs but got %d', ...
              s, numel(valid), nFound);
    end
    
    % 3) build a full Ns×K matrix, stuffing missing ones with NaN
    full_means = nan(nSubj, nROIs);
    full_means(:, valid) = roi_means;
    
    % 4) store
    parc_data{s} = full_means;
end

