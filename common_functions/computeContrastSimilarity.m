function simMat = computeContrastSimilarity(contrastCell, labels)
% COMPUTECONTRASTSIMILARITY  Correlate Nc study‐contrast patterns and plot
%
%   simMat = computeContrastSimilarity(contrastCell, labels)
%
% INPUTS
%   contrastCell – 1×Nc cell, each cell is [Ns×Nv] (subjects×voxels)
%   labels       – (opt) 1×Nc cell array of strings for each study
%
% OUTPUT
%   simMat       – Nc×Nc correlation matrix of mean‐pattern similarity
%
% USAGE:
%   simMat = computeContrastSimilarity(contrastCell, studyNames);

Nc = numel(contrastCell);
if nargin<2, labels = arrayfun(@num2str, 1:Nc, 'Uni',false); end

% 1) build Nc×Nv matrix of mean patterns
Nv = size(contrastCell{1},2);
meanPatterns = nan(Nc, Nv);
for i = 1:Nc
    M = contrastCell{i};        % [Ns×Nv]
    meanPatterns(i,:) = mean(M, 1, 'omitnan');
end

% 2) compute pairwise Pearson corr
simMat = corr(meanPatterns');   % → Nc×Nc

% 3) plot heatmap
figure;
imagesc(simMat, [-1 1]);
colormap(flipud(hot)); colorbar;
set(gca, ...
    'XTick', 1:Nc, 'XTickLabel', labels, ...
    'YTick', 1:Nc, 'YTickLabel', labels);
xtickangle(45);
title('Study‐wise pattern similarity (r)','Interpreter','none');

% 4) (optional) hierarchical clustering dendrogram
Z = linkage(1-simMat, 'average');
figure;
dendrogram(Z, 0, 'Labels', labels, 'Orientation','left');
title('Cluster tree of studies by pattern dissimilarity','Interpreter','none');

end
