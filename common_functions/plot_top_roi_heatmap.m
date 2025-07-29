function plot_top_roi_heatmap(study_means, region_labels,  study_names,threshold, topN)
% PLOT_TOP_ROI_HEATMAP  Rank ROIs by prevalence and plot a heatmap
%
% INPUTS:
%   study_means   – nStudies×nROIs matrix (from average_subjects)
%   region_labels – 1×nROIs cell of ROI names (atlas.labels)
%   threshold     – scalar, t‑score threshold for “active” 
%   topN          – number of top ROIs to display
%
% This will:
%  1. Count for each ROI how many studies exceed threshold.
%  2. Sort ROIs by that count (descending).
%  3. Plot an imagesc heatmap of the topN ROIs.

% 1) compute “active‑study” counts per ROI
counts = sum(abs(study_means) > threshold, 1);

% 2) rank and select
[~, sortIdx] = sort(counts, 'descend');
topIdx = sortIdx(1:min(topN, numel(sortIdx)));

% 3) extract data & labels
top_data   = study_means(:, topIdx);
top_labels = region_labels(topIdx);

% 4) plot
figure;
imagesc(top_data);
colorbar;
xlabel('ROI');
ylabel('Study');
set(gca, ...
    'XTick', 1:numel(topIdx), ...
    'XTickLabel', top_labels, ...
    'XTickLabelRotation', 90, ...
    'YTick', 1:size(study_means,1), ...
    'YTickLabel',  study_names, ...
    'TickLabelInterpreter','none');

title(sprintf('Top %d ROIs by studies > %.2f t‑score', numel(topIdx), threshold));
axis tight;
end
