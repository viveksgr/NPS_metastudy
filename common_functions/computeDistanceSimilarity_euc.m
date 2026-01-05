function simMat = computeDistanceSimilarity_euc(contrastCell, labels)
% COMPUTEDISTANCESIMILARITY  Second‐order similarity via ROI–ROI distances
%
%   simMat = computeDistanceSimilarity(contrastCell, labels)
%
% INPUTS
%   contrastCell – 1×Nc cell, each cell is [Ns×Nv] (subjects×ROIs)
%   labels       – (optional) 1×Nc cell of study names
%
% OUTPUT
%   simMat       – Nc×Nc Spearman correlation of vectorized ROI–ROI distance mats

Nc = numel(contrastCell);
if nargin<2 || isempty(labels)
    labels = arrayfun(@(i) sprintf('Study%d',i), 1:Nc, 'Uni',false);
end

% Pre‐allocate
Nv = size(contrastCell{1},2);
ut   = triu(true(Nv),1);
vlen = nnz(ut);
distVecs = nan(Nc, vlen);

for i = 1:Nc
    % 1) average across subjects → 1×Nv
    M = contrastCell{i};
    mu = mean(M, 1, 'omitnan');     % 1×Nv

    % 2) compute Nv×Nv Euclidean distance matrix
    %    pdist returns 1×(Nv*(Nv-1)/2) vector of distances already in ut order
    dvec = pdist(mu', 'euclidean'); % 1×vlen

    distVecs(i,:) = dvec;
end

% 3) compute Nc×Nc Spearman correlation across those distance‐vectors
simMat = corr(distVecs', distVecs', 'Type','Spearman');

% 4) plot heatmap
figure;
% a = logical(eye(size(simMat)));
% simMat(a) = nan;
imagesc(simMat, [-1 1]);
% colormap(flipud(hot)); 
colorbar;
caxis([0 0.4])
set(gca, ...
    'XTick', 1:Nc, 'XTickLabel', labels, ...
    'YTick', 1:Nc, 'YTickLabel', labels, ...
    'TickLabelInterpreter','none');
xtickangle(45);
title('Study–Study Similarity (Spearman of ROI–ROI distances)', ...
      'Interpreter','none');
end
