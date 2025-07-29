mean_map = average_subjects_tcorrected(parc_data_multi);

% --- 2) Build a plain image_vector from ds.volInfo + your mean_map
iv = image_vector( ...
       'volInfo', ds.volInfo, ...
       'dat',     mean_map     );









montage(ms_brain{1}, ...
    'threshold',   2.5, ...        % e.g. show t > 2.5
    'montagetype','full', ...      % 'compact','multirow', etc.
    'color',      'hot', ...
    'noblobs','nooutline');        % turn off extra overlays

dat = ms_brain{1};
o = fmridisplay(ms_brain{1}, 'display_mode','montage');
o = addblobs(o, dat>3, 'color',[0 0 1]);   % overlay suprathreshold mask


% make sure SPM defaults are loaded:
spm('Defaults','FMRI');

% then
orthviews('Reset');
orthviews(dat);


%% Plotting
% --- 1) Suppose “ds” is your fmri_data object with
%         ds.dat = [nVoxels × nSubjects]

ds = ms_brain{1};
% average across subjects:
mean_map = mean(ds.dat, 2, 'omitnan');   % [nVoxels×1] vector

% --- 2) Build a plain image_vector from ds.volInfo + your mean_map
iv = image_vector( ...
       'volInfo', ds.volInfo, ...
       'dat',     mean_map     );

% --- 3) (Re)load SPM defaults so montage won’t crash on tol_orient
spm('Defaults','FMRI');

% --- 4) Show a slice montage of your group-mean image
montage(iv, ...
    'threshold',   0, ...        % show all values > 0
    'montagetype','full', ...    % axial/cor/sag slices
    'color',       'hot', ...    % colormap
    'noblobs', ...               % turn off ROI outlines/blobs
    'nooutline');                % (redundant with noblobs)

% or, if you prefer the fmridisplay engine directly:
h = fmridisplay(iv, 'display_mode','montage');
colorbar;
title('Mean activation across subjects');



%%
% 1) Make sure SPM defaults (tol_orient etc.) are in memory
spm('Defaults','FMRI');

% 2) Clear out any baseline/surface info on your image_vector
iv.surface = [];    % so montage won’t try render_on_surface
iv.baseimgs = {};   % so montage won’t try to read a T1

% 3) Call montage with everything turned off except the overlay itself
montage(iv, ...
    'threshold',   0, ...       % show everything > 0
    'montagetype','full', ...   % axial/cor/sag slices
    'color',       'hot', ...   % your colormap
    'nobaseline', ...           % no anatomical underlay
    'noblobs', ...              % no ROI‐outlines or extra blobs
    'nooutline');             

vol = nan(iv.volInfo.dim);
vol(iv.volInfo.wh_inmask) = iv.dat;
niftiwrite(vol, 'roi_map.nii', iv.volInfo);


% % 1) Suppose you already have:
% %    atlas        – your canlab atlas object
% %    roi_vals     – 1×NROI vector of values you want to plot
% 
% % 2) Make a voxel-wise map
% codes      = atlas.dat;                  % V×1 integer ROI codes
% voxel_map  = zeros(size(codes));         % initialize V×1 map
% for r = 1:numel(atlas.labels)
%     voxel_map(codes==r) = roi_vals(r);
% end
% 
% % 3) Package into an image_vector *manually*
% iv = image_vector( ...
%        'volInfo', atlas.volInfo, ...     % header/dimensions/wh_inmask
%        'dat',     voxel_map );           % your V×1 data vector
% 
% spm('Defaults','FMRI');     % or 'EEG' / whatever modality you’re using
% spm_jobman('initcfg');      % if you haven’t already
global defaults
defaults.images.tol_orient = 1e-4;
% 
% montage(iv, ...
%     'threshold',   0, ...         % show everything > 0
%     'montagetype','full', ...
%     'color',       'hot', ...
%     'noblobs', ...                % <-- turn off all blob overlays
%     'nooutline');                % <-- supresses atlas borders if any