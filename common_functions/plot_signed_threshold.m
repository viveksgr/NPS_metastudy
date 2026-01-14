function o2 = plot_signed_threshold(tiv_rs, T, varargin)
% PLOT_SIGNED_THRESHOLD_COLORMAP  show ±T map with graded diverging colormap
%
% o2 = plot_signed_threshold_colormap(tiv_rs, T, 'transvalue',0.25, 'ncolors',256)
%
% - tiv_rs: image_vector (V×1) of t or effect values
% - T:      positive threshold (e.g. tinv)
% Optional name/value:
%   'transvalue' (default 0.2)  transparency for overlay
%   'ncolors'    (default 256)  total colormap size (will split across neg/pos)
%
% Returns: canlab_results_fmridisplay object (o2)

% parse optional args
p = inputParser;
addParameter(p,'transvalue',0.8,@isnumeric);
addParameter(p,'ncolors',256,@(x) isnumeric(x) && isscalar(x) && x>=8);
parse(p,varargin{:});
tv = p.Results.transvalue;
ncolors = p.Results.ncolors;

% basic checks
if ~isa(tiv_rs,'image_vector')
    error('tiv_rs must be a canlab image_vector');
end

dat = double(tiv_rs.dat(:));   % Vx1 ensure double
dat(isnan(dat)) = 0;

% mask null band
masked = dat;
masked(abs(masked) < T) = 0;

if all(masked==0)
    warning('No voxels exceed ±%g. Nothing to show.', T);
    o2 = [];
    return
end

% symmetric range for diverging map
maxAbs = max(abs(masked(:)));
cmapHalf = round(ncolors/2);

% define colour-stops (you can tweak these RGB triplets)
darkBlue = [0.00 0.06 0.40];   % deep blue
cyan     = [0.00 0.70 0.85];
midGrey  = [0.80 0.80 0.80];   % neutral centre (near zero)
orange   = [0.92 0.45 0.07];
yellow   = [1.00 0.88 0.14];

% make neg colormap: darkBlue -> cyan -> midGrey
negStops = [darkBlue; cyan; midGrey];
tneg = linspace(0,1,cmapHalf)';
cneg = interp1([0;0.5;1], negStops, tneg, 'linear');

% make pos colormap: midGrey -> orange -> yellow
posStops = [midGrey; orange; yellow];
tpos = linspace(0,1,ncolors-cmapHalf)';
cpos = interp1([0;0.5;1], posStops, tpos, 'linear');

% combine to one continuous colormap (neg first, pos second)
cmap = [cneg; cpos];

% create new image_vector with masked values (same volInfo)
iv = image_vector('volInfo', tiv_rs.volInfo, 'dat', masked);

% create display object
o2 = canlab_results_fmridisplay([], 'nobaseline', 'nooutline', 'full');

% Try to add as a single diverging layer (preferred for surface rendering)
try
    % note: 'split' tells the renderer to treat negative and positive sides separately
    o2 = addblobs(o2, region(iv, 'noverbose'), ...
                  'split', ...
                  'cmap', cmap, ...
                  'cmaprange', [-maxAbs maxAbs], ...
                  'transvalue', tv, ...
                  'nooutline');
    % build legend
    o2 = legend(o2, 'noverbose');
    return
catch ME
    % if the single-pass approach fails (different versions of canlab/addblobs), fall back
    warning('Single-layer diverging colormap failed (%s). Falling back to two-layer rendering.', ME.message);
end

% FALLBACK: separate negative and positive layers (graded but separate)
posdat = masked; posdat(posdat < 0) = 0;
negdat = masked; negdat(negdat > 0) = 0;

% pos colormap uses second half, neg colormap uses first half
if any(negdat(:) ~= 0)
    rngneg = double([min(negdat(negdat~=0)), -T]); % min (neg) to threshold
    try
        o2 = addblobs(o2, region(image_vector('volInfo', tiv_rs.volInfo, 'dat', negdat), 'noverbose'), ...
                      'mincolor', darkBlue, ...
                      'cmap', cneg, ...
                      'cmaprange', [min(negdat(negdat~=0)) -T], ...
                      'transvalue', tv, 'nooutline');
    catch
        o2 = addblobs(o2, region(image_vector('volInfo', tiv_rs.volInfo, 'dat', negdat), 'noverbose'), ...
                      'mincolor', darkBlue, 'transvalue', tv, 'nooutline');
    end
end

if any(posdat(:) ~= 0)
    rngpos = double([T, max(posdat(posdat~=0))]);
    try
        o2 = addblobs(o2, region(image_vector('volInfo', tiv_rs.volInfo, 'dat', posdat), 'noverbose'), ...
                      'maxcolor', yellow, ...
                      'cmap', cpos, ...
                      'cmaprange', [T max(posdat(posdat~=0))], ...
                      'transvalue', tv, 'nooutline');
    catch
        o2 = addblobs(o2, region(image_vector('volInfo', tiv_rs.volInfo, 'dat', posdat), 'noverbose'), ...
                      'maxcolor', yellow, 'transvalue', tv, 'nooutline');
    end
end

o2 = legend(o2, 'noverbose');
end
