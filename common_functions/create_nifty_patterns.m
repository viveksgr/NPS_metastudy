function create_nifty_patterns(atlas,mean_data,fname)
%% Brain plot
spm('Defaults','FMRI');  

% loads SPM’s global defaults
% atlas = load_atlas('canlab2023_coarse_fmriprep20_1mm');
% → atlas.dat       is [V×1] integer ROI codes (1…N)
% → atlas.labels    is {N×1} cell array of names
% → atlas.volInfo   is an SPM-style header (dim, mat, etc.)

N = numel(atlas.labels);          % number of parcels
if numel(mean_data) ~= N
    error('mean_data length (%d) ≠ atlas parcels (%d).', numel(mean_data), N);
end

voxel_vec = zeros(size(atlas.dat));      % V×1 output vector
for r = 1:N
    voxel_vec(atlas.dat == r) = mean_data(r);
end

vol3d = zeros(atlas.volInfo.dim);         % 3-D array in atlas space
vol3d(atlas.volInfo.wh_inmask) = voxel_vec;

V            = atlas.volInfo;    % copy header
V.fname      = fullfile(pwd, fname);   % output filename
V.descrip    = 'Mean meta-analysis value per ROI';
% Tell SPM to save the data as 32-bit float instead of the atlas’ int type
V.dt         = [spm_type('float32') 0];       % 16 → float32
V.pinfo      = [1; 0; 0];                     % no scaling
spm_write_vol(V, vol3d);         % writes uncompressed .nii

% templ_gz = fullfile( ...
%     fileparts(which('canlab_config')), ...
%     'canlab_canonical_brains', 'Canonical_brains_surfaces', ...
%     'fmriprep20_template.nii.gz');
% 
% gunzip(templ_gz);                      % creates fmriprep20_template.nii
% 
% 
% spm_orthviews('Reset');
% T = spm_vol('fmriprep20_template.nii');
% M = spm_vol('roi_meta_map.nii');
% spm_orthviews('Image', T);
% spm_orthviews('AddBlobs',1, M, ones(size(vol3d)), '');   % overlay
% spm_orthviews('Redraw');