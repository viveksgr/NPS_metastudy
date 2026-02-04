function [fig, order, allMeans] = plotInterventionROIHeat(dataCell, varargin)
% plotInterventionROIHeat  Heatmap of ROI×Intervention means, interventions clustered
%
% USAGE
%   [fig, order, allMeans] = plotInterventionROIHeat(dataCell, ...)
%
% INPUTS
%   dataCell : 1 x Ncell cell array. Each cell is (nSubjects x nROIs). Subjects may vary per cell.
%
% Name-Value options:
%   'InterventionLabels' (default {}) : cellstr labels for interventions (length N)
%   'Distance'           (default 'correlation') : distance metric for pdist (used on intervention profiles)
%   'Linkage'            (default 'average') : linkage method for hierarchical clustering
%   'Reorder'            (default true) : whether to reorder interventions by clustering
%   'Scale'              (default 'none') : 'none' | 'zscoreROIs' (zscore across interventions per ROI)
%   'Clim'               (default []) : color limits [min max]; if empty force symmetric around 0
%   'FigureTitle'        (default '') : title string
%
% OUTPUTS
%   fig      - figure handle
%   order    - permutation of interventions (1..N)
%   allMeans - nROIs x nInterventions matrix of means used for plotting
%
% Example:
%   [fig,order,allMeans] = plotInterventionROIHeat(dataCell, 'InterventionLabels', labels);
%

% parse inputs
p = inputParser;
addParameter(p,'InterventionLabels',{}, @(x) iscellstr(x) || isstring(x) || isempty(x));
addParameter(p,'Distance','correlation', @ischar);
addParameter(p,'Linkage','average', @ischar);
addParameter(p,'Reorder', true, @islogical);
addParameter(p,'Scale','none', @(s) ismember(s, {'none','zscoreROIs'}));
addParameter(p,'Clim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(p,'FigureTitle','', @ischar);
parse(p,varargin{:});
labels = cellfun(@char, cellstr(p.Results.InterventionLabels), 'UniformOutput', false);
distanceMetric = p.Results.Distance;
linkageMethod  = p.Results.Linkage;
doReorder      = p.Results.Reorder;
scaleMode      = p.Results.Scale;
clim_in        = p.Results.Clim;
figTitle       = p.Results.FigureTitle;

% basic checks
N = numel(dataCell);
if N == 0, error('dataCell is empty'); end

% compute means across subjects for each cell (use nanmean to handle NaNs)
for i = 1:N
    Xi = dataCell{i};
    if ~ismatrix(Xi)
        error('Each cell must be a 2D matrix (subjects x ROIs). Problem at cell %d', i);
    end
    % compute mean across rows = subjects -> 1 x nROIs
    m = nanmean(double(Xi), 1);
    if i == 1
        nROIs = numel(m);
    else
        if numel(m) ~= nROIs
            error('All cells must have the same number of ROIs (columns). Cell %d mismatch', i);
        end
    end
    allMeans(:, i) = m(:);   %#ok<AGROW> nROIs x N
end

% optional scaling by ROI (zscore across interventions)
switch lower(scaleMode)
    case 'none'
        % keep as-is
    case 'zscorerois'
        allMeans = (allMeans - mean(allMeans,2))./ (std(allMeans,0,2) + eps);
    otherwise
        error('Unknown scale mode %s', scaleMode);
end

% compute ordering of interventions by clustering (on columns)
order = 1:N;
if doReorder && N > 1
    try
        Dvec = pdist(allMeans','correlation'); % distance between intervention profiles
        Z = linkage(Dvec, linkageMethod);
        % try optimal leaf ordering first
        try
            order = optimalleaforder(Z, Dvec);
        catch
            % fallback: build a hidden dendrogram and get permutation
            h = figure('Visible','off');
            [H,T,order] = dendrogram(Z, 0); %#ok<ASGLU>
            close(h);
        end
    % catch ME
    %     warning('Clustering failed (%s). Using original order.', ME.message);
    %     order = 1:N;
    end
end

% prepare labels
if isempty(labels)
    labels = arrayfun(@(k) sprintf('Int%d', k), 1:N, 'Uni', false);
end
labels = labels(:);
if numel(labels) ~= N
    warning('Number of labels does not match number of interventions. Using generic labels.');
    labels = arrayfun(@(k) sprintf('Int%d', k), 1:N, 'Uni', false);
end

% colormap: blue-white-red
ncol = 256;
half = floor(ncol/2);
c1 = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];       % blue->white
c2 = [ones(ncol-half,1), linspace(1,0,ncol-half)', linspace(1,0,ncol-half)']; % white->red
cmap = [c1; c2];

% color limits
if isempty(clim_in)
    mx = max(abs(allMeans(:)));
    clim = [-mx mx];
else
    clim = clim_in;
end

% Plot heatmap
fig = figure('Color','w', 'Name','ROI x Intervention heatmap');
imagesc(1:N, 1:nROIs, allMeans(:, order)', clim);
colormap(cmap);
colorbar;
set(gca, 'YDir', 'normal'); % ensure ROI1 at top or bottom preference
xlabel('Intervention'); ylabel('ROI');
title(figTitle, 'Interpreter','none');

% ticks
set(gca, 'XTick', 1:N, 'XTickLabel', labels(order), 'XTickLabelRotation', 45, 'TickLabelInterpreter','none');
set(gca, 'YTick', 1:min(20,nROIs)); % show subset of ROI ticks if too many
% if many ROIs, better to not label all; user can tweak

% nice plotting defaults
axis tight;
box on;

% output
allMeans = allMeans; %# already nROIs x N, return for convenience

end
