function [mu_meta, var_meta] = average_subjects_tcorrected(parc_data)
% parc_data: cell array of size [1 x nStudies]
% each cell is [N_s x nROIs] of subject‐wise contrasts

nStudies = numel(parc_data);
nROIs    = size(parc_data{1},2);

% Preallocate accumulators
sum_wy = zeros(1, nROIs);   % sum of w_i * mean_i
sum_w  = zeros(1, nROIs);   % sum of w_i

epsVar = 1e-12;  % small threshold to avoid dividing by zero

for s = 1:nStudies
    data = parc_data{s};          % [N_s x nROIs]
    Ns   = size(data,1);

    % Compute per‐ROI mean and variance
    mu_s   = nanmean(data,1);     % 1 x nROIs
    var_s  = nanvar(data,0,1);    % 1 x nROIs (sample variance)

    % Variance of the mean
    var_mean = var_s ./ Ns;       % 1 x nROIs

    % Avoid zero or near-zero variances
    var_mean(var_mean < epsVar) = Inf;

    % Inverse‐variance weights (zero for zero‐variance ROIs)
    w = 1 ./ var_mean;            % 1 x nROIs

    % Accumulate for meta‐analysis
    sum_wy = sum_wy + (w .* mu_s);
    sum_w  = sum_w  + w;
end

% Final meta‐analytic estimates per ROI
mu_meta  = sum_wy ./ sum_w;     % 1 x nROIs
var_meta = 1   ./ sum_w;        % 1 x nROIs  (inverse of sum of weights)
end
