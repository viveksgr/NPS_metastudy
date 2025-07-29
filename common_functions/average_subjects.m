function study_means = average_subjects(parc_data)
% AVERAGE_SUBJECTS  Collapse each study’s subject × ROI matrix to a study × ROI mean
%
% INPUT:
%   parc_data – 1×nStudies cell, each cell is [nSubjects × nROIs] of t‑scores
%
% OUTPUT:
%   study_means – nStudies×nROIs matrix of mean t‑score per ROI (ignoring NaNs)

nStudies = numel(parc_data);
nROIs    = size(parc_data{1},2);
study_means = nan(nStudies, nROIs);

for s = 1:nStudies
    % mean across subjects, omit NaNs for missing ROIs
    study_means(s,:) = mean(parc_data{s}, 1, 'omitnan');
end
end
