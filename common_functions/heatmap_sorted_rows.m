function [fig, order, clusterIdx] = heatmap_sorted_rows(M, varargin)
% HEATMAP_SORTED_ROWS  Heatmap of nStudy x nROI matrix with rows sorted by hierarchical clustering
%
%   [fig, order, clusterIdx] = heatmap_sorted_rows(M)
%   [fig, order, clusterIdx] = heatmap_sorted_rows(M, 'RowLabels', rlab, 'ColLabels', clav)
%
% INPUTS
%   M          - numeric matrix (nStudy x nROI)
% Name-value pairs (optional):
%   'RowLabels' - cellstr or string array, length nStudy (labels for rows after sorting)
%   'ColLabels' - cellstr or string array, length nROI  (labels for columns)
%
% OUTPUTS
%   fig        - figure handle
%   order      - permutation of row indices (1..nStudy) after sorting (so M(order,:) is plotted)
%   clusterIdx - vector length nStudy with values 1 or 2 assigning each study to cluster 1 or 2
%
% Notes:
% - Uses pdist(...,'correlation') and linkage(...,'average').
% - If optimalleaforder is available it's used; otherwise falls back to dendrogram ordering.

% ----------------- parse inputs -----------------------
p = inputParser;
addParameter(p,'RowLabels',{}, @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(p,'ColLabels',{}, @(x) isempty(x) || iscellstr(x) || isstring(x));
parse(p, varargin{:});
rowLabels = p.Results.RowLabels;
colLabels = p.Results.ColLabels;

[nStudy, nROI] = size(M);
if nStudy < 2
    error('Need at least 2 rows (studies) to cluster.');
end

% ---------------- clustering rows ---------------------
% distance between rows (treat rows as feature vectors across ROIs)
D = pdist(M, 'correlation');        % 1 x (nchoosek(nStudy,2))
Z = linkage(D, 'average');          % hierarchical tree

% try to get optimal leaf order; fallback to dendrogram order (quietly)
try
    order = optimalleaforder(Z, D);
catch
    % dendrogram returns order but creates a figure; suppress by making invisible
    h = figure('Visible','off');
    [~,~,order] = dendrogram(Z, 0);
    close(h);
end

% ---------------- two-cluster assignment ----------------
% cluster into 2 groups
clusterIdx = cluster(Z, 'maxclust', 2); % vector length nStudy, values in {1,2}
% Optionally reorder the labels so cluster 1 corresponds to top half of dendrogram:
% We keep clusterIdx as returned.

% ---------------- plot heatmap -------------------------
fig = figure('Color','w');
imagesc(M(order, :)); 
axis tight;
colormap(parula);
colorbar;
xlabel('ROIs');
ylabel('Studies (clustered)');

% set x / y ticks / labels if provided
if ~isempty(colLabels)
    if numel(colLabels) ~= nROI
        warning('ColLabels length mismatch; ignoring col labels.');
    else
        set(gca, 'XTick', 1:nROI, 'XTickLabel', colLabels, 'XTickLabelRotation', 45, 'TickLabelInterpreter','none');
    end
else
    set(gca, 'XTick', []); % avoid crowding if many ROIs
end

if ~isempty(rowLabels)
    if numel(rowLabels) ~= nStudy
        warning('RowLabels length mismatch; ignoring row labels.');
    else
        % show labels in clustered order
        rl = rowLabels(order);
        set(gca, 'YTick', 1:nStudy, 'YTickLabel', rl, 'TickLabelInterpreter','none');
    end
end

title('Heatmap (rows clustered)');
set(gca, 'YDir', 'normal'); % put first entry at top

% ---------------- return clusterIdx aligned with original rows -----------
% clusterIdx is already aligned to original row order; if user wants cluster
% indices in the plotted order, they can use clusterIdx(order)
end
