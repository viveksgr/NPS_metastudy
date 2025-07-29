function group_t = average_subjects_t(parc_data)
% Runs a one‐sample t‐test (vs. zero) for each study × ROI
%
% INPUT:
%   parc_data : 1×nStudies cell array, where each cell is [N_s × nROIs]
%
% OUTPUT:
%   group_t   : nStudies × nROIs matrix of t‐statistics

nStudies = numel(parc_data);
nROIs    = size(parc_data{1},2);

% Preallocate
group_t = nan(nStudies, nROIs);

for s = 1:nStudies
    data = parc_data{s};          % [N_s × nROIs]
    for r = 1:nROIs
        y = data(:,r);
        y = y(~isnan(y));         % drop NaNs
        if numel(y) > 1
            [~,~,~,stats] = ttest(y, 0);
            group_t(s,r) = stats.tstat;
        end
    end
end
end
