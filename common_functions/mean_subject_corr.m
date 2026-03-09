function [mean_r, se_r, mean_z, se_z, nPairs] = mean_subject_corr(datobj, varargin)
% MEAN_SUBJECT_CORR  Mean pairwise correlation among subjects (and SE)
%
% Usage:
%   [mean_r, se_r] = mean_subject_corr(datobj)
%   [mean_r, se_r, mean_z, se_z, nPairs] = mean_subject_corr(datobj, 'useFisher', true, 'verbose', true)
%
% Inputs:
%   datobj        - canlab object with .dat (nVox x nSubjects) OR numeric matrix (nVox x nSubjects)
%
% Name/value options:
%   'useFisher'   - (true) compute mean & SE in Fisher z space (recommended)
%   'verbose'     - (false)
%
% Outputs:
%   mean_r  - mean Pearson r across all subject pairs (upper triangle mean)
%   se_r    - approximate SE of mean_r (z-space SE back-propagated)
%   mean_z  - mean Fisher z (atanh of r) used for computations
%   se_z    - SE in z-space: std(z)/sqrt(nPairs)
%   nPairs  - number of subject pairs used (nSubs * (nSubs-1) / 2 minus NaN pairs)
%
% Notes:
%  - Uses pairwise-complete correlations: voxels with NaNs are ignored for each subject pair.
%  - For stability, correlations are clamped to (-0.999999, 0.999999) before atanh.

% parse options
p = inputParser;
addParameter(p,'useFisher', true, @islogical);
addParameter(p,'verbose', false, @islogical);
parse(p, varargin{:});
useFisher = p.Results.useFisher;
doVerbose = p.Results.verbose;

% get data matrix X: nV x nS
if isstruct(datobj) || isobject(datobj)
    if isprop(datobj,'dat') || isfield(datobj, 'dat')
        X = datobj.dat;
    else
        error('Provided object does not have a .dat field.');
    end
else
    X = datobj;
end

if ~isnumeric(X)
    error('Input data must be numeric or object with numeric .dat field.');
end

% ensure double
X = double(X);
[nV, nS] = size(X);
if nS < 2
    warning('Need at least 2 subjects to compute pairwise correlations.');
    mean_r = NaN; se_r = NaN; mean_z = NaN; se_z = NaN; nPairs = 0;
    return;
end

% compute subject×subject Pearson correlation matrix using pairwise-complete
% columns are subjects, so corr(X, 'Rows','pairwise') computes corr among columns
R = corr(X, 'Rows', 'pairwise');

% extract upper triangle (exclude diagonal)
upperMask = triu(true(nS), 1);
rvals = R(upperMask);

% handle NaNs that may appear if a pair has <2 non-NaN voxels
valid = ~isnan(rvals);
rvals = rvals(valid);
nPairs = numel(rvals);

if nPairs == 0
    warning('No valid subject pairs (all correlations are NaN).');
    mean_r = NaN; se_r = NaN; mean_z = NaN; se_z = NaN;
    return;
end

% clamp extreme values to avoid infinite atanh
epsClamp = 1 - 1e-6;
rvals(rvals >=  epsClamp) =  epsClamp;
rvals(rvals <= -epsClamp) = -epsClamp;

% Fisher z transform
zvals = atanh(rvals);   % Fisher z

% compute mean & SE
mean_z = mean(zvals, 'omitnan');
se_z   = std(zvals, 0, 'omitnan') ./ sqrt(nPairs);

if useFisher
    % back-transform mean z to r, and propagate SE approx via derivative
    mean_r = tanh(mean_z);
    % approximate se in r-space using delta method: dr/dz = 1 - r^2
    se_r = (1 - mean_r.^2) * se_z;
else
    % compute mean directly in r-space and naive SE
    mean_r = mean(rvals, 'omitnan');
    se_r   = std(rvals, 0, 'omitnan') ./ sqrt(nPairs);
end

if doVerbose
    fprintf('Subjects: %d, pairs used: %d\n', nS, nPairs);
    fprintf('Mean r = %.4f, SE_r = %.4f (mean_z = %.4f, SE_z = %.4f)\n', mean_r, se_r, mean_z, se_z);
end

end
