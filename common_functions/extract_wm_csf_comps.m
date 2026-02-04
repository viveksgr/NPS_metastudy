function wm_csf_comps = extract_wm_csf_comps(dat,sw)
if nargin<2
    sw = true;
end
%
[value, components] = extract_gray_white_csf(dat, 'masks', ...
    {'gray_matter_mask.nii', 'canonical_white_matter_thrp5_ero1.nii', ...
    'canonical_ventricles_thrp5_ero1.nii'});
wm_nuisance = value(:,2:3); % white, CSF
wm_nuisance_comps = [components{2} components{3}]; % 5 components for white, CSF
if sw
    wm_csf_comps = [wm_nuisance scale(wm_nuisance_comps)];
else
    wm_csf_comps = [wm_nuisance];
end
