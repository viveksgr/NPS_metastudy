function plotMDSFromSimilarity(simMat, C, labels)
% PLOTMDSFROMSIMILARITY  MDS plot of studies from a similarity matrix
%
%   plotMDSFromSimilarity(simMat, C, labels)
%
% INPUTS:
%   simMat – Nc×Nc similarity matrix (e.g. Spearman r between studies)
%   C      – Nc×1 binary or two‐level numeric vector (e.g. 0/1 or 1/2)
%             Determines dot color: values==1 → red, values==2 → blue
%   labels – (opt) Nc×1 cell array of study names
%
% This function performs classical MDS on the distance matrix D = 1 - simMat,
% extracts the first two dimensions, and scatter‐plots the studies in 2D,
% coloring each point red or blue depending on C.

Nc = size(simMat,1);
if nargin<3 || isempty(labels)
    labels = arrayfun(@(i) sprintf('S%d',i), 1:Nc, 'Uni',false);
end
if numel(C)~=Nc
    error('Color code vector C must have length %d', Nc);
end

% 1) Convert similarity to distance
D = 1 - simMat;
% ensure zero diagonal
D(1:Nc+1:end) = 0;

% 2) Classical MDS (cmdscale)
[Y, eigvals] = cmdscale(D);
% Y is Nc×Nc, but we only need first 2 dims
coords = Y(:,1:2);

% 3) Plot
figure; hold on;
% define colors: assume C has values 1 and 2
clr = zeros(Nc,3);
clr(C==1,:) = repmat([1 0 0], sum(C==1), 1);  % red
clr(C==2,:) = repmat([0 0 1], sum(C==2), 1);  % blue

scatter(coords(:,1), coords(:,2), 80, clr, 'filled');
text(coords(:,1), coords(:,2), labels, ...
     'VerticalAlignment','bottom', 'HorizontalAlignment','right','Interpreter','none');

xlabel('MDS Dimension 1');
ylabel('MDS Dimension 2');
title('MDS of Study Similarity','Interpreter','none');
grid on; axis equal;
hold off;
end
