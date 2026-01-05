function simMat = computeContrastSimilarity_rsa(contrastCell, labels)
% COMPUTECONTRASTSIMILARITY  Study–study similarity via PCA + Spearman
%
%   simMat = computeContrastSimilarity(contrastCell, labels)
%
% INPUTS
%   contrastCell – 1×Nc cell, each cell is [Ns×Nv] (subjects×voxels)
%   labels       – (optional) 1×Nc cell array of study names
%
% OUTPUT
%   simMat       – Nc×Nc matrix of Spearman correlations between studies
%
% This version does a “second‐order” reduction by PCA (retain components
% explaining ≥90% variance, up to max 10 PCs), then computes Spearman
% correlation between the reduced study profiles.

Nc = numel(contrastCell);
if nargin<2 || isempty(labels)
    labels = arrayfun(@(i) sprintf('Study%d',i), 1:Nc, 'Uni',false);
end

% 1) build Nc×V matrix of mean patterns
V = size(contrastCell{1},2);
M = nan(Nc, V);
for i = 1:Nc
    M(i,:) = mean(contrastCell{i}, 1, 'omitnan');
end

% 2) PCA reduction
[~, score, ~, ~, explained] = pca(M, 'Rows','complete');
cumVar = cumsum(explained);
K = find(cumVar >= 80, 1, 'first');    % components to retain ≥90% var
% K = min(K, 10);                        % cap at 10 PCs
R = score(:, 1:K);                     % Nc×K reduced matrix

% 3) compute Spearman correlation between rows of R
simMat = corr(R', R', 'Type', 'Spearman');

% 4) plot heatmap
figure('Position',[0.5 0.5 720 640]);
imagesc(simMat, [-1 1]);
% colormap("jet");
colorbar;
set(gca, ...
    'XTick', 1:Nc, 'XTickLabel', labels, ...
    'YTick', 1:Nc, 'YTickLabel', labels, ...
    'TickLabelInterpreter','none');
xtickangle(90);
title('Study–Study Similarity (Spearman on PCA scores)', 'Interpreter','none');
end
