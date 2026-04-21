function report = diagnose_Y_matrix(Y, varargin)
% DIAGNOSE_Y_MATRIX  Diagnostic checks on the N_subjects x N_voxels data matrix Y.
%
%   report = diagnose_Y_matrix(Y)
%   report = diagnose_Y_matrix(Y, 'nvox_sample', 5000, 'verbose', true)
%
% Checks for numerical pathologies that cause discretization or chaos in OLS:
%   1. Data type and precision
%   2. Integer-valued data (int16 nifti source quantization)
%   3. Repeating scale factor (suggests integer data with slope scaling)
%   4. Effective unique values per voxel (low = quantized)
%   5. Per-voxel variance distribution (zero-variance columns)
%   6. NaN / Inf count
%   7. Condition number of Y'*Y (collinearity across subjects)
%   8. Subject-wise mean distribution (inter-subject heterogeneity)
%   9. Histogram of a random voxel sample (visual check)
%
% INPUTS
%   Y               - N_subjects x N_voxels numeric matrix
%
% OPTIONAL NAME/VALUE
%   'nvox_sample'   - voxels to subsample for expensive checks (default 5000)
%   'verbose'       - print summary to command window (default true)
%
% OUTPUT
%   report          - struct with all diagnostic results

p = inputParser;
addRequired(p, 'Y', @isnumeric);
addParameter(p, 'nvox_sample', 5000, @isscalar);
addParameter(p, 'verbose', true, @islogical);
parse(p, Y, varargin{:});

nv_samp  = p.Results.nvox_sample;
verbose  = p.Results.verbose;

[N, V] = size(Y);
report = struct();
report.N = N;
report.V = V;

sep = repmat('-', 1, 60);

if verbose, fprintf('\n%s\nDIAGNOSE_Y_MATRIX  [%d subjects x %d voxels]\n%s\n', sep, N, V, sep); end

% -------------------------------------------------------------------------
% 1. Data type
% -------------------------------------------------------------------------
report.class = class(Y);
if verbose, fprintf('1. Data type        : %s\n', report.class); end
if ~isa(Y, 'double')
    if verbose, fprintf('   WARNING: Y is not double. OLS should be run in double precision.\n'); end
end

% -------------------------------------------------------------------------
% 2. NaN / Inf
% -------------------------------------------------------------------------
report.n_nan = sum(isnan(Y(:)));
report.n_inf = sum(isinf(Y(:)));
if verbose
    fprintf('2. NaN count        : %d  (%.2f%%)\n', report.n_nan, 100*report.n_nan/numel(Y));
    fprintf('   Inf count        : %d  (%.2f%%)\n', report.n_inf, 100*report.n_inf/numel(Y));
end

% -------------------------------------------------------------------------
% 3. Integer-valued check (int16 source quantization)
%    If most values are whole numbers, data came from int16 nifti
% -------------------------------------------------------------------------
Yf = double(Y(:));
Yf = Yf(isfinite(Yf));
frac_integer = mean(Yf == round(Yf));
report.frac_integer_valued = frac_integer;
if verbose
    fprintf('3. Integer-valued   : %.1f%% of finite values\n', 100*frac_integer);
    if frac_integer > 0.5
        fprintf('   WARNING: >50%% integer values — likely int16 nifti source.\n');
        fprintf('   Quantization will discretize OLS outputs.\n');
    end
end

% -------------------------------------------------------------------------
% 4. Scale factor fingerprint
%    For int16 data with slope: values = int16 * slope + intercept
%    The minimum nonzero absolute difference reveals the scale step
% -------------------------------------------------------------------------
vox_idx = randperm(V, min(nv_samp, V));
Ysamp = double(Y(:, vox_idx));
diffs = diff(sort(Ysamp(:)));
diffs = diffs(diffs > 0 & isfinite(diffs));
if ~isempty(diffs)
    min_step = min(diffs);
    report.min_nonzero_step = min_step;
    if verbose
        fprintf('4. Min nonzero step : %g  (if ~integer, may be int16 scale factor)\n', min_step);
    end
else
    report.min_nonzero_step = NaN;
end

% -------------------------------------------------------------------------
% 5. Effective unique values per voxel (subsample)
%    Low counts per voxel = quantization
% -------------------------------------------------------------------------
n_unique = zeros(1, numel(vox_idx));
for k = 1:numel(vox_idx)
    col = Ysamp(:, k);
    n_unique(k) = numel(unique(col(isfinite(col))));
end
report.median_unique_per_vox = median(n_unique);
report.min_unique_per_vox    = min(n_unique);
if verbose
    fprintf('5. Unique vals/vox  : median=%g  min=%g  (max possible=%d subjects)\n', ...
        report.median_unique_per_vox, report.min_unique_per_vox, N);
    if report.median_unique_per_vox < N * 0.5
        fprintf('   WARNING: Few unique values per voxel — strong quantization signal.\n');
    end
end

% -------------------------------------------------------------------------
% 6. Per-voxel variance distribution
% -------------------------------------------------------------------------
col_var = var(double(Y), 0, 1);  % 1 x V
report.n_zero_var_vox  = sum(col_var == 0);
report.n_low_var_vox   = sum(col_var < 1e-6);
report.median_col_var  = median(col_var);
if verbose
    fprintf('6. Zero-var voxels  : %d  |  Near-zero-var (<1e-6): %d\n', ...
        report.n_zero_var_vox, report.n_low_var_vox);
    fprintf('   Median vox var   : %g\n', report.median_col_var);
end

% -------------------------------------------------------------------------
% 7. Per-subject mean and variance (detect outlier subjects / global shifts)
% -------------------------------------------------------------------------
row_mean = mean(double(Y), 2);   % N x 1
row_var  = var(double(Y), 0, 2); % N x 1
report.subject_means = row_mean;
report.subject_vars  = row_var;
[~, outlier_mean_idx] = sort(abs(row_mean - median(row_mean)), 'descend');
if verbose
    fprintf('7. Subject means    : median=%.3g  MAD=%.3g\n', ...
        median(row_mean), mad(row_mean, 1));
    fprintf('   Most deviant subjects (by mean): rows %s\n', ...
        num2str(outlier_mean_idx(1:min(5,N))'));
end

% -------------------------------------------------------------------------
% 8. Condition of Y (via singular values of subsample)
%    High condition -> collinear subjects -> unstable OLS
% -------------------------------------------------------------------------
Ysamp_d = double(Ysamp);
Ysamp_d(~isfinite(Ysamp_d)) = 0;
sv = svd(Ysamp_d, 'econ');
sv = sv(sv > 0);
report.condition_number = sv(1) / sv(end);
if verbose
    fprintf('8. Condition (Y sub): %.3g  (>1e6 suggests near-collinear subjects)\n', ...
        report.condition_number);
end

% -------------------------------------------------------------------------
% 9. Histogram plots
% -------------------------------------------------------------------------
fig = figure('Position', [100 100 900 320], 'Name', 'diagnose_Y_matrix');

% Panel 1: distribution of all sampled values
subplot(1,3,1);
histogram(Ysamp(:), 100, 'EdgeColor', 'none');
title('Value distribution (vox sample)');
xlabel('Value'); ylabel('Count');

% Panel 2: per-voxel variance (log scale)
subplot(1,3,2);
histogram(log10(col_var(col_var > 0)), 80, 'EdgeColor', 'none');
title('log_{10}(per-voxel var)');
xlabel('log_{10}(variance)'); ylabel('Count');

% Panel 3: per-subject mean
subplot(1,3,3);
plot(sort(row_mean), 'o', 'MarkerSize', 4);
title('Per-subject global mean (sorted)');
xlabel('Subject rank'); ylabel('Mean');

sgtitle('diagnose\_Y\_matrix');
report.fig = fig;

if verbose, fprintf('%s\n', sep); end

end
