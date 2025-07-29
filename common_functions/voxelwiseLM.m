function [tmap_iv, pmap_iv, df] = voxelwiseLM(data_cell, C, nuisance_cell)
% VOXELWISELM   Voxel‐wise linear model across studies + nuisance
%
%   [tmap_iv, pmap_iv, df] = voxelwiseLM(data_cell, C, nuisance_cell)
%
% INPUTS
%   data_cell    – 1×S cell, each fmri_data.dat = [V×n_s]
%   C            – 1×S contrast weights for each study
%   nuisance_cell– (opt) 1×S cell, each [n_s×p] covariate matrix
%
% OUTPUTS
%   tmap_iv      – image_vector of voxel‐wise t‐stats for the study contrast
%   pmap_iv      – image_vector of two‐tailed p‐values
%   df           – degrees of freedom = N - rank(X)

if nargin<3, nuisance_cell = cell(size(data_cell)); end
S = numel(data_cell);
assert(numel(C)==S,'Contrast length must match number of studies');

% 1) Stack data & build regressors
V = data_cell{1}.volInfo.n_inmask;
Y = [];        % will be [N×V]
G = [];        % study contrast regressor [N×1]
Z = [];        % nuisances [N×p]
for s=1:S
    D = data_cell{s}.dat;        % [V×n_s]
    ns = size(D,2);
    Y  = [Y; D'];                % append [n_s×V]
    G  = [G; C(s)*ones(ns,1)];
    Zs = nuisance_cell{s};
    if isempty(Zs)
        Z = [Z; zeros(ns,0)];  
    else
        assert(size(Zs,1)==ns,'nuisance rows must match subjects');
        Z = [Z; Zs];
    end
end
[N, ~] = size(Y);

% 2) Design matrix X = [1 G Z]
X = [ones(N,1), G, Z];
df = N - rank(X);

% 3) OLS betas for each voxel
XtXinv = inv(X'*X);               % (p+2)×(p+2)
B      = XtXinv * (X' * Y);       % (p+2)×V

% 4) Residual MSE per voxel
E   = Y - X*B;                     
mse = sum(E.^2,1) / df;           % 1×V

% 5) t‐stat for the study‐contrast regressor (column 2 of B)
varB2 = mse * XtXinv(2,2);        % 1×V
tvals = B(2,:) ./ sqrt(varB2);

% 6) two‐tailed p‐values
pvals = 2*tcdf(-abs(tvals), df);
[~,p_val2] = fdr_benjhoc(pvals);

% 7) wrap as image_vectors
volInfo = data_cell{1}.volInfo;
tmap_iv = image_vector('volInfo',volInfo,'dat',tvals(:));
pmap_iv = image_vector('volInfo',volInfo,'dat',pvals(:));

montage(tmap_iv, 'trans',        0.5,'montagetype','full');

tiv = tmap_iv;
T = tinv(1-0.025, df);

% display as a full slice montage, red blobs, semi-transparent, no baseline
montage(tmap_iv, ...
    'threshold',    T, ...            % only |t|>=T
    'cmaprange',    [T max(tmap_iv.dat)], ...  % color only positive t’s
    'color',        [1 .2 .2], ...    % red colormap
    'trans',        0.5, ...          % 50% transparency
    'montagetype', 'full', ...
    'nobaseline', ...
    'noblobs',      false, ...        % allow blobs
    'nooutline');                      % drop ROI edges
end


% 1) load the underlay as an image_vector
tmpl = image_vector('image_names','fmriprep20_template.nii');

% 2) resample your t-map onto it
tiv_rs = resample_space(tiv, tmpl);

% 3) now montage with no worries of misalignment
montage(tiv_rs, ...
    'threshold',T, 'trans',0.3, ...
    'montagetype','full');