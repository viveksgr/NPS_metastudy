function resampled_list = alignFmriDataToReference(fmri_list, refIdx)
% ALIGNFMRIDATATOREFERENCE  Resample a list of fmri_data objects to one common grid
%
%   resampled_list = alignFmriDataToReference(fmri_list, refIdx)
%
% INPUTS:
%   fmri_list – 1×N cell array of fmri_data objects
%   refIdx    – which index in fmri_list to use as the reference grid
%                (must be integer between 1 and N)
%
% OUTPUT:
%   resampled_list – 1×N cell array of fmri_data objects,
%                    each now on the same volInfo grid as fmri_list{refIdx}
%
% EXAMPLE:
%   % align all studies to the first one
%   aligned = alignFmriDataToReference(fmriStudies, 1);

N = numel(fmri_list);
if refIdx < 1 || refIdx > N
    error('refIdx must be between 1 and %d', N);
end

% pick out the reference
refObj = fmri_list{refIdx};

% build a “template” image_vector that carries only the volInfo of refObj
templateIV = image_vector( ...
    'volInfo', refObj.volInfo, ...
    'dat',     zeros(refObj.volInfo.n_inmask,1) );

resampled_list = cell(1, N);

for i = 1:N
    src = fmri_list{i};
    switch compare_space(src, templateIV)
      case 0
        % identical space & mask → no change
        resampled_list{i} = src;

      case 1
        % same physical space but different sampling → resample!
        resampled_list{i} = resample_space(src, templateIV);

      otherwise
        error('Study %d is in an incompatible space; cannot resample.', i);
    end
end
end
