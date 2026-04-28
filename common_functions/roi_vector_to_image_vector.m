function iv = roi_vector_to_image_vector(atlas, roi_vals)
% ROI_VECTOR_TO_IMAGE_VECTOR  Map per-ROI values back into voxel space
%
%   iv = roi_vector_to_image_vector(atlas, roi_vals)
%
%   atlas    - canlab atlas object (e.g. load_atlas('canlab2024_coarse_...'))
%   roi_vals - numeric vector, length == atlas.n_regions, one value per ROI
%              (NaN or 0 for regions you want left blank)
%
%   Returns an image_vector with the same volInfo as the atlas, suitable
%   for passing directly to plot_signed_threshold().

if ~isa(atlas, 'atlas')
    error('atlas must be a canlab atlas object');
end
if numel(roi_vals) ~= atlas.num_regions
    error('roi_vals must have %d elements (one per atlas region), got %d', ...
          atlas.num_regions, numel(roi_vals));
end

labels = double(atlas.dat(:));   % V×1 integer label per voxel (0 = background)
vox_vals = zeros(size(labels));  % default to 0 (transparent in plot)

for r = 1:atlas.num_regions
    vox_vals(labels == r) = roi_vals(r);
end

% NaN → 0 so plot_signed_threshold masks them cleanly
vox_vals(isnan(vox_vals)) = 0;

iv = image_vector('volInfo', atlas.volInfo, 'dat', vox_vals);
end
