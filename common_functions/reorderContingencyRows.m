function OUT2 = reorderContingencyRows(OUT, varargin)
% reorderContingencyRows  Reorder study rows of OUT.contingencyProps by hierarchical clustering
%
% OUT2 = reorderContingencyRows(OUT)
% OUT2 = reorderContingencyRows(OUT, 'K',2, 'Distance','correlation', 'Linkage','average', 'Plot',true)
%
% Inputs:
%  - OUT: struct returned from studyClusterCrosstab (must contain .contingencyProps)
%
% Optional name/value pairs:
%  - 'K'         : number of clusters to cut the dendrogram into (default 2)
%  - 'Distance'  : pdist distance metric (default 'correlation')
%  - 'Linkage'   : linkage method (default 'average')
%  - 'Plot'      : true/false show dendrogram + heatmap (default true)
%  - 'LeafOrder' : provide an external leaf ordering (vector of indices) to use instead of recomputing
%
% Output OUT2 contains the same fields as OUT, plus:
%  - OUT2.reorder_index         : permutation of studies applied (row indices)
%  - OUT2.contingencyProps_re  : reordered contingencyProps
%  - OUT2.contingencyCounts_re : reordered contingencyCounts
%  - OUT2.studyIDs_re          : reordered studyIDs
%  - OUT2.clusterAssignments   : cluster assignment per study (after reordering)
%  - OUT2.linkage              : linkage matrix used
%  - OUT2.leafOrder            : leaf order used to reorder rows
%
% Notes:
%  - This function clusters rows (studies) by their contingency proportion vectors.
%  - If you pass 'LeafOrder', that order is used directly (bypasses recompute).

% --- parse inputs
p = inputParser;
addParameter(p,'K',2,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'Distance','correlation',@(x)ischar(x)||isstring(x));
addParameter(p,'Linkage','average',@(x)ischar(x)||isstring(x));
addParameter(p,'Plot',true,@islogical);
addParameter(p,'LeafOrder',[], @(x)isnumeric(x) || isempty(x));
parse(p,varargin{:});
K = p.Results.K;
distMetric = char(p.Results.Distance);
linkMethod = char(p.Results.Linkage);
doplot = p.Results.Plot;
leafOrderUser = p.Results.LeafOrder;

% basic checks
if ~isfield(OUT,'contingencyProps')
    error('OUT must have field .contingencyProps created by studyClusterCrosstab.');
end

props = OUT.contingencyProps;
counts = OUT.contingencyCounts;
studyIDs = OUT.studyIDs;

[nStudies, nClusters] = size(props);

if ~isempty(leafOrderUser)
    if numel(leafOrderUser) ~= nStudies
        error('Provided LeafOrder length must equal number of studies (%d).', nStudies);
    end
    leafOrder = leafOrderUser(:);
    Z = [];
else
    % compute distance and linkage
    if nStudies <= 1
        leafOrder = 1:nStudies;
        Z = [];
    else
        % if the chosen distance metric is correlation, pdist returns 1-corr by default
        D = pdist(props, distMetric);
        Z = linkage(D, linkMethod);
        % try to find a leaf order that puts similar rows nearby
        try
            leafOrder = optimalleaforder(Z, D);
        catch
            % fallback: use dendrogram to produce leaf order (creates a figure momentarily)
            try
                h = figure('Visible','off'); 
                [~,~,leafOrder] = dendrogram(Z, 0);
                close(h);
            catch
                % last fallback: use the naive order from linkage (1..N)
                leafOrder = 1:nStudies;
            end
        end
    end
end

% reorder outputs
props_re = props(leafOrder, :);
counts_re = counts(leafOrder, :);
studyIDs_re = studyIDs(leafOrder);

% cluster assignments on the (reordered) linkage (if we have a linkage)
if ~isempty(leafOrderUser) && isempty(Z)
    % compute Z from data to get cluster assignment
    D = pdist(props, distMetric);
    Z = linkage(D, linkMethod);
end

if ~isempty(Z)
    clusters_all = cluster(Z, 'maxclust', K);
else
    clusters_all = ones(nStudies,1);
end
% reorder cluster labels according to leafOrder
clusters_re = clusters_all(leafOrder);

% assemble OUT2 (copy original then add/replace fields)
OUT2 = OUT;
OUT2.reorder_index = leafOrder;
OUT2.contingencyProps_re = props_re;
OUT2.contingencyCounts_re = counts_re;
OUT2.studyIDs_re = studyIDs_re;
OUT2.clusterAssignments = clusters_re;
OUT2.linkage = Z;
OUT2.leafOrder = leafOrder;

% if original assignedGroup/binaryGroup exist, reorder those too
if isfield(OUT,'assignedGroup')
    OUT2.assignedGroup_re = OUT.assignedGroup(leafOrder);
end
if isfield(OUT,'binaryGroup')
    OUT2.binaryGroup_re = OUT.binaryGroup(leafOrder);
end

% optional plotting: dendrogram (left) + heatmap (right)
if doplot
    hf = figure('Color','w','Name','Study x Cluster proportions (reordered)');
    % layout: left narrow dendro, right heatmap
    t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
    % dendrogram
    nexttile(1,[1 1]);
    if ~isempty(Z)
        % show dendrogram; reorder by leafOrder for consistent appearance
        [H, T, ~] = dendrogram(Z, 0, 'Reorder', leafOrder);
        title('Dendrogram (studies)');
        ylabel('Distance');
        set(gca,'XTick',[]);
    else
        text(0.5,0.5,'No dendrogram (single study)','HorizontalAlignment','center');
        axis off
    end
    % heatmap
    nexttile(2,[1 1]);
    imagesc(props_re);
    colormap(gca, parula);
    colorbar;
    xlabel('Clusters');
    ylabel('Studies (reordered)');
    title('Per-study proportion into clusters');
    set(gca,'YTick',1:nStudies,'YTickLabel',studyIDs_re,'TickLabelInterpreter','none');
    % if many studies, reduce tick labels
    if nStudies > 60
        set(gca,'YTick',1:floor(nStudies/10):nStudies);
    end
end

end
