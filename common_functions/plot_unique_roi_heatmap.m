function plot_unique_roi_heatmap( group_t, region_labels,  study_names)

% assume:
%   group_t       = [nStudies × K]  matrix of study-level t-scores
%   region_labels = 1×K cell array of ROI names (atlas.labels)
%
thresh = 3.1;
nStudies = size(group_t,1);

% 1) find ROIs exceeding threshold in exactly one study
binmat   = group_t > thresh;                  % logical [S×K]
uniqROI  = find(sum(binmat,1)==1);            % indices of ROIs unique to one study

% 2) for each unique ROI, record which study and its t-value
nUniq    = numel(uniqROI);
stud     = zeros(1,nUniq);
tvals    = zeros(1,nUniq);
for k = 1:nUniq
    r = uniqROI(k);
    i = find(binmat(:,r),1,'first');          % the one study
    stud(k)  = i;
    tvals(k)  = group_t(i,r);
end

% 3) sort those ROIs by descending t-value
[~,ord]  = sort(tvals,'descend');
uniqROI  = uniqROI(ord);
stud     = stud(ord);

% 4) pick top 20 (or fewer if less than 20 unique)
N = min(20, numel(uniqROI));
selROI = uniqROI(1:N);

% 5) build heatmap data
H = group_t(:, selROI);  % [nStudies×N], each column hot in exactly one row

% 6) plot
figure;
imagesc(H);
colorbar;
set(gca, ...
    'XTick', 1:N, ...
    'XTickLabel', region_labels(selROI), ...
    'XTickLabelRotation', 45, ...
    'YTick', 1:nStudies, ...
   'YTickLabel',  study_names, ...
    'TickLabelInterpreter','none' ...
);
title(sprintf('Top %d ROIs uniquely > %.1f t-score', N, thresh));
axis tight;
