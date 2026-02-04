function save_image_vector_as_nifti(iv, out_fname)
% iv : image_vector or fmri_data object
% out_fname : full path, e.g. '/tmp/mymap.nii'

% reconstruct_image returns a 3D (or 4D) volume for the in-mask voxels
[voldata, vectorized_voldata, xyz_struct] = reconstruct_image(iv);
% voldata: X x Y x Z x Images  (if multiple images)
vol3d = voldata(:, :, :, 1);   % choose first image (or loop images)

% find a template header: try iv.volInfo.V if available, otherwise fall back
if isfield(iv, 'volInfo') && isfield(iv.volInfo, 'V') && ~isempty(iv.volInfo.V)
    V = iv.volInfo.V;   % V may be a struct or cell; use first if array
    if iscell(V), V = V{1}; end
else
    % fallback: use a canonical template shipped with CanlabCore (adjust path)
    V = spm_vol(which('SPM8_colin27T1_seg.img')); % or any NIfTI template on path
end

% prepare header for output
V.fname = out_fname;
V.dt = [16 0];              % float32
V.dim = size(vol3d);
% spm_write_vol expects the data in same orientation/shape as V
spm_write_vol(V, vol3d);
end
