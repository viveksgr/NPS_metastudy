function [tmap_iv, pmap_iv, df,t_thr] = voxelwiseLM(data_cell, C, st_vec, regnum,nuisance_cell)
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
% regnum = 1;
if nargin<5, nuisance_cell = cell(size(data_cell)); end
S = numel(data_cell);
assert(numel(C)==S,'Contrast length must match number of studies');

% 1) Stack data & build regressors
V = data_cell{1}.volInfo.n_inmask;
Y = [];        % will be [N×V]
G = [];        % study contrast regressor [N×1]
Z = [];        % nuisances [N×p]
M = [];       % Contrast specific reg
M2 = [];        % Study specific reg
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

    % Add study_spec regressors
    Ms = ones(ns,1);
    M = blkdiag(M, Ms);

    if s==1
        M2 = Ms;
    else
        cdr = st_vec(s)==st_vec(s-1);
        if cdr
            nrs = size(M2,2)-1;
            M2 = [M2; [zeros(ns,nrs) Ms]];
        else
            M2 = blkdiag(M2, Ms);
        end
    end
end
[N, ~] = size(Y);

% Column reduction:
% idxs = find(diff(diff(st_vec))==-1)+1;
% M2 = M2(:,st_vec(idxs));

% 2) Design matrix X = [1 G Z]
X = [ones(N,1), G, Z, M]; % Change to add custom regressors
% X = [ones(N,1), Z, M];

% Orthog the study design
% [X, ~] = makeInterceptGrandMean(X, 3:size(X,2));
studyCols = 2:size(X,2); %3:size(X,2); % Change to add custom regressors

[X, ~] = makeInterceptGrandMean(X, studyCols, 'keepIntercept', true, 'frontRegs', false);

% X = [ones(N,1), G];
df = N - rank(X);

% 3) OLS betas for each voxel
% XtXinv = inv(X'*X);               % (p+2)×(p+2)
% B      = XtXinv * (X' * Y);       % (p+2)×V
B      = double(pinv(X) * Y);  

% 4) Residual MSE per voxel
E      = Y - X*B;                     
mse    = sum(E.^2,1) / df;           % 1×V

% 5) Estimate covariance of beta‐hat via pinv(X'X)
XtX_pinv = pinv(X'*X);               % (p×p) pseudoinverse

% 6) t‐stat for the study‐contrast regressor (column 2 of B)
%    (or whichever column index your contrast lives in)
col       =regnum;                       % e.g. column 2 for G regressor
varB_col  = mse .* XtX_pinv(col,col); % 1×V
tvals     = B(col,:) ./ sqrt(varB_col);

% 6) two‐tailed p‐values
pvals = 2*tcdf(-abs(tvals), df);
[p1,~] = fdr_benjhoc(pvals);

t_thr(1) = tinv(1-fdr_benjhoc( tcdf((tvals), df)),df);
t_thr(2) = tinv(1-fdr_benjhoc( tcdf(-(tvals), df)),df);

% 7) wrap as image_vectors

volInfo = data_cell{1}.volInfo;
tmap_iv = image_vector('volInfo',volInfo,'dat',tvals(:));
pmap_iv = image_vector('volInfo',volInfo,'dat',pvals(:));

tmpl = image_vector('image_names','fmriprep20_template.nii');
tiv_rs = resample_space(tmap_iv, tmpl);
% montage(tiv_rs, 'trans', 0.5,'montagetype','full');

% T = tinv(1-0.025, df);
T = tinv(1-p1,df);
% 
% save_image_vector_as_nifti(tiv_rs, 'statmap.nii')     
% % save stat image to nifti
% 
% % run MatlabTFCE (example; follow that repo's API)
% tfce_img = matlab_tfce('statmap.nii');
% 
% % load TFCE result into canlab
% tfce_iv = image_vector('image_names','tfce_statmap.nii');



plot_signed_threshold(tiv_rs, T)
% montage(tiv_rs, 'trans', 0.2,'montagetype','full','threshold',    T,'cmaprange',    [min(tmap_iv.dat) max(tmap_iv.dat)]);
% 
% % display as a full slice montage, red blobs, semi-transparent, no baseline
% montage(tiv_rs, ...
%     'threshold',    T, ...            % only |t|>=T
%     'cmaprange',    [T max(tmap_iv.dat)], ...  % color only positive t’s
%     'color',        [1 .2 .2], ...    % red colormap
%     'trans',        0.5, ...          % 50% transparency
%     'montagetype', 'full', ...
%     'nobaseline', ...
%     'noblobs',      false, ...        % allow blobs
%     'nooutline');                      % drop ROI edges
end

% st_vec = [1 2 3 4 4 4 5 6 6 6 6 7 7 8 9 10];
% C = [1 1 -1 -1 -1 -1 1 1 1 1 -1 -1 -1 -1 1];
% % 1) load the underlay as an image_vector
% tmpl = image_vector('image_names','fmriprep20_template.nii');
% 
% % 2) resample your t-map onto it
% tiv_rs = resample_space(tiv, tmpl);
% 
% % 3) now montage with no worries of misalignment
% montage(tiv_rs, ...
%     'threshold',T, 'trans',0.3, ...
%     'montagetype','full');