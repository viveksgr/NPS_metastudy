function [data_clean, keptIdx, removedIdx, stats] = remove_low_variance_features(data, varargin)
% REMOVE_LOW_VARIANCE_FEATURES Remove nearly-constant columns from data.
%   data: observations x features
%   Optional name/value:
%     'RelSTDThreshold' (default 1e-3) - remove if std(col)/median(std(cols)) < this
%     'AbsSTDThreshold' (default eps)  - also remove if std(col) < this
%     'Names' (default {})             - cellstr of feature names (1 x nfeatures)
%     'Verbose' (default true)
%
% Returns cleaned data, kept and removed indices, and stats struct.

p = inputParser;
addParameter(p,'RelSTDThreshold',1e-3,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'AbsSTDThreshold',eps,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'Names',{},@iscell);
addParameter(p,'Verbose',true,@islogical);
parse(p,varargin{:});
rt = p.Results.RelSTDThreshold;
at = p.Results.AbsSTDThreshold;
names = p.Results.Names;
verbose = p.Results.Verbose;

% compute std across rows for each column
colstd = std(data, 0, 1, 'omitnan');   % 1 x nfeatures
medstd = median(colstd(colstd>0));     % median among nonzero stds to be robust
if isempty(medstd) || medstd==0
    medstd = mean(colstd) + eps;
end
relstd = colstd ./ (medstd + eps);

% choose removal mask
removeMask = (relstd < rt) | (colstd < at);
keptIdx = find(~removeMask);
removedIdx = find(removeMask);

data_clean = data(:, keptIdx);

stats.colstd = colstd;
stats.relstd = relstd;
stats.medstd = medstd;
stats.removed_count = numel(removedIdx);

if verbose
    fprintf('remove_low_variance_features: removed %d / %d features (relstd < %g or std < %g)\n', ...
        numel(removedIdx), numel(colstd), rt, at);
    if ~isempty(names) && ~isempty(removedIdx)
        fprintf('Removed feature names (first 10):\n');
        disp(names(removedIdx(1:min(10,numel(removedIdx)))));
    end
end
end
