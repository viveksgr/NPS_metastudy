function [pc_scores, study_means2D] = plot_study_kmeans(parc_data, k, group_idx)
% PLOT_STUDY_KMEANS  Cluster subjects across multiple studies and plot 2-D study means
%
%   [pc_scores, study_means2D] = plot_study_kmeans(parc_data, k, group_idx)
%
% INPUTS
%   parc_data   1×S cell, each cell = [nSubj_s × nROI] t-score matrix
%   k           desired number of k-means clusters
%   group_idx   1×S numeric labels assigning each study to a user-defined group
%               (e.g. [1 1 1 2 2 1 1 2 2 1]).  Pass [] or omit for default coloring.
%
% OUTPUTS
%   pc_scores       [Ntotal×2] PCA space for every subject (for custom plots)
%   study_means2D   [S×2] PCA coordinates of each study mean
%
% EXAMPLE
%   plot_study_kmeans(parc_data, 4, [1 1 2 2 1 1 2 2 1 1 2]);
%
%   — Vivek & ChatGPT, 2025-06-12

% ----------------------------------------------------------
% 0.  Handle inputs, defaults
% ----------------------------------------------------------
if nargin < 3 || isempty(group_idx)
    group_idx = 1:numel(parc_data);    % unique color per study
end
S = numel(parc_data);                  % number of studies
if numel(group_idx) ~= S
    error('group_idx must have same length (%d) as parc_data.', S);
end

% ----------------------------------------------------------
% 1.  Concatenate all subjects
% ----------------------------------------------------------
allData   = [];
study_tag = [];                        % study index per subject
for s = 1:S
    allData   = [allData; parc_data{s}];          %#ok<AGROW>
    study_tag = [study_tag; s * ones(size(parc_data{s},1),1)]; %#ok<AGROW>
end
[Ntotal, nROI] = size(allData);

% ----------------------------------------------------------
% 2.  z-score features across subjects, run k-means
% ----------------------------------------------------------
% Z = zscore(allData);                   % N×nROI standardized

% ----------------------------------------------------------
% 2.  PCA → dimensionality reduction, then k-means
% ----------------------------------------------------------
allData(isnan(allData))=0;
Z = zscore(allData);                       % N×nROI



% ---- PCA ----
[coeff, score, latent, ~, explained] = pca(Z,'Rows','pairwise', ...  % keep partial rows
                                           'Algorithm','svd');      % robust for tall matrices

% choose how many PCs to keep
cumexp  = cumsum(explained);
nPC     = min( find(cumexp >= 90, 1, 'first'), 10 );  % ≥90 % var OR ≤10 PCs
reduced = score(:, 1:nPC);               % N×nPC  (data for clustering)

% ---- k-means on reduced space ----
[idx, ~] = kmeans(reduced, k, 'Replicates',20, 'Display','off');

% ----------------------------------------------------------
% 3.  Keep first two PCs for plotting
% ----------------------------------------------------------
pc_scores = score(:, 1:2);               % N×2  (for scatter)


% compute study means in that 2-D space
study_means2D = nan(S,2);
for s = 1:S
    study_means2D(s,:) = mean(pc_scores(study_tag==s, :), 1, 'omitnan');
end

% ----------------------------------------------------------
% 4.  Plot
% ----------------------------------------------------------
figure; hold on;
colors = lines(max(group_idx));        % distinct color per user group

% (a) light subject points
scatter(pc_scores(:,1), pc_scores(:,2), 12, .8*[1 1 1], 'filled', ...
        'MarkerEdgeAlpha',0.1, 'MarkerFaceAlpha',0.2);

% (b) study means in stronger color
for s = 1:S
    g = group_idx(s);
    plot(study_means2D(s,1), study_means2D(s,2), ...
        'o', 'MarkerSize',10, 'MarkerFaceColor',colors(g,:), ...
        'MarkerEdgeColor','k', 'LineWidth',1.1);
    text(study_means2D(s,1), study_means2D(s,2), sprintf('  S%02d',s), ...
        'FontSize',8, 'Color','k', 'VerticalAlignment','middle');
end

% (c) cluster centroids (for reference)
cent = zeros(k,2);
for c = 1:k
    cent(c,:) = log(abs(mean(pc_scores(idx==c,:), 1, 'omitnan')));
end
plot(cent(:,1), cent(:,2), 'kx', 'MarkerSize',12, 'LineWidth',2);

xlabel('PC-1'); ylabel('PC-2');
grid on; axis square;
title(sprintf('k-means clustering (%d clusters) in ROI space', k));
legend({'subjects','study means','k-means centroids'}, 'Location','bestoutside');
end
