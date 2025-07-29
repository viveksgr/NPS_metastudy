function [tmap_iv, pmap_iv, df] = fixedEffectsContrast(data_cell, C, alpha)
% FIXEDEFFECTSCONTRAST  Voxel‐wise fixed‐effects contrast across studies
%
%   [tmap_iv, pmap_iv, df] = fixedEffectsContrast(data_cell, C, alpha)
%
% INPUTS
%   data_cell : 1×S cell of fmri_data objects, each dat=[V×n_s]
%   C          : 1×S contrast weights (sum does *not* have to be zero)
%   alpha      : (opt) two‐tailed significance level, default = 0.05
%
% OUTPUTS
%   tmap_iv    : image_vector of length-V t‐statistic map
%   pmap_iv    : image_vector of length-V two‐tailed p‐value map
%   df         : degrees of freedom (uses min(n_s−1) across studies)
%
% Example:
%   C = [1 -1];
%   [tiv, piv, df] = fixedEffectsContrast({d1,d2}, C);
%   montage(tiv, 'threshold', tinv(1-0.025,df));

if nargin<3, alpha = 0.05; end
S = numel(data_cell);
if numel(C)~=S
    error('Contrast length (%d) must match number of studies (%d).', numel(C), S);
end

% 1) Gather per‐study stats
V = data_cell{1}.volInfo.n_inmask;  % number of voxels in the mask
means    = zeros(S, V);
vars     = zeros(S, V);
ns       = zeros(1, S);
for s = 1:S
    dat = data_cell{s}.dat;         % [V × n_s]
    ns(s) = size(dat,2);
    means(s,:) = mean(dat, 2, 'omitnan')';
    vars(s,:)  = var( dat, 0, 2, 'omitnan')';
end

% 2) Compute contrast numerator & variance
num = C * means;    % 1×V
den = sqrt( sum( (C(:).^2) .* (vars ./ ns(:)), 1 ) );  % 1×V

% avoid divide-by-zero
den(den==0) = eps;
tvals = num ./ den;  % 1×V

% 3) Degrees of freedom (conservative)
df = min(ns)-1;

% 4) Two‐tailed p‐values
pvals = 2 * tcdf(-abs(tvals), df);
[~,p_val2] = fdr_benjhoc(pvals);

% 5) Wrap into image_vector
volInfo = data_cell{1}.volInfo;    % common header
tmap_iv  = image_vector('volInfo', volInfo, 'dat', tvals(:));
pmap_iv  = image_vector('volInfo', volInfo, 'dat', pvals(:));

% 6) (Optional) threshold display example
montage(tmap_iv, 'threshold', tinv(1-alpha/2,df), 'montagetype','full');

end
