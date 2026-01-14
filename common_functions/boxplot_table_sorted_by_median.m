function h = boxplot_table_sorted_by_median(T, group, varargin)
%BOXPLOT_TABLE_SORTED_BY_MEDIAN Boxplot each table column (NaN-padded),
%sorted by increasing median, with same-group boxes in same color.
%
%   h = boxplot_table_sorted_by_median(T, group)
%   h = boxplot_table_sorted_by_median(T, group, 'Name', Value, ...)
%
% Inputs
%   T      : table (each variable/column is one box; NaN padding allowed)
%   group  : grouping variable with length = width(T) (numeric/cellstr/string/categorical)
%
% Name-Value options
%   'Ax'          : target axes handle (default = gca)
%   'Palette'     : Kx3 colormap for K groups (default = lines(K))
%   'BoxWidth'    : scalar passed to boxplot (default = 0.6)
%   'ShowOutliers': true/false (default = true)
%   'LabelMode'   : 'vars' or 'none' (default = 'vars')  % x tick labels
%
% Output
%   h : struct of handles returned by boxplot (plus some extras)

% -------- parse inputs
p = inputParser;
p.addRequired('T', @(x) istable(x));
p.addRequired('group', @(x) numel(x) == width(T));
p.addParameter('Ax', [], @(x) isempty(x) || isgraphics(x,'axes'));
p.addParameter('Palette', [], @(x) isempty(x) || (ismatrix(x) && size(x,2)==3));
p.addParameter('BoxWidth', 0.6, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('ShowOutliers', true, @(x) islogical(x) && isscalar(x));
p.addParameter('LabelMode', 'vars', @(x) any(strcmpi(x,{'vars','none'})));
p.addParameter('SortMode', 'column', ...
    @(x) any(strcmpi(x,{'column','group','none'})));

p.parse(T, group, varargin{:});

ax = p.Results.Ax;
if isempty(ax), ax = gca; end

% -------- extract data as numeric matrix (nRows x nCols)
X = table2array(T);
if ~isnumeric(X)
    error('T must contain numeric variables only.');
end
nCols = size(X,2);

% -------- normalize group to categorical (for stable mapping)
g = group;
if isstring(g) || iscellstr(g)
    g = string(g);
end
g = categorical(g);

% -------- compute medians per column ignoring NaNs
med = nan(1,nCols);
for j = 1:nCols
    med(j) = median(X(~isnan(X(:,j)), j), 'omitnan');
end

% Columns with all-NaN => median becomes NaN; push them to the end
sortMode = lower(p.Results.SortMode);

switch sortMode
    case 'none'
        ord = 1:nCols;

    case 'column'
        medSortKey = med;
        medSortKey(isnan(medSortKey)) = +Inf;
        [~, ord] = sort(medSortKey, 'ascend');

    case 'group'
        % Compute median for each group
        cats0 = categories(g);
        G = numel(cats0);
        groupMed = nan(G,1);

        for k = 1:G
            idx = g == cats0{k};
            groupMed(k) = median(med(idx), 'omitnan');
        end

        % Sort groups by their median
        [~, gOrd] = sort(groupMed, 'ascend');
        sortedCats = cats0(gOrd);

        % Within each group, preserve column order (or optionally sort by med)
        ord = [];
        for k = 1:G
            idx = find(g == sortedCats{k});
            ord = [ord, idx(:)']; %#ok<AGROW>
        end
end


Xo  = X(:,ord);
go  = g(ord);




% -------- palette per group
cats = categories(go);
K = numel(cats);
pal = p.Results.Palette;
if isempty(pal)
    pal = lines(K);
else
    if size(pal,1) < K
        error('Provided Palette has %d rows but needs at least %d (one per group).', size(pal,1), K);
    end
    pal = pal(1:K,:);
end

% -------- build long-form vectors for boxplot
xPos = repelem(1:nCols, size(Xo,1));         % positions, repeated down rows
yVal = Xo(:);

% remove NaNs so boxplot doesn't see them as data points
keep = ~isnan(yVal);
xPos = xPos(keep);
yVal = yVal(keep);

% -------- draw boxplot
cla(ax);
hold(ax,'on');

sym = 'o';
if ~p.Results.ShowOutliers
    sym = '';  % boxplot: empty symbol => no outliers plotted
end

hBox = boxplot(ax, yVal, xPos, ...
    'Positions', 1:nCols, ...
    'Widths', p.Results.BoxWidth, ...
    'Symbol', sym);

% Labels
switch lower(p.Results.LabelMode)
    case 'vars'
        varNames = string(T.Properties.VariableNames);
        varNames = varNames(ord);
        ax.XTick = 1:nCols;
        ax.XTickLabel = varNames;
        ax.XTickLabelRotation = 45;
    case 'none'
        ax.XTick = 1:nCols;
        ax.XTickLabel = [];
end

ax.XLim = [0.5, nCols+0.5];
grid(ax,'on');

% -------- color each box according to group
% boxplot returns boxes as line objects with Tag 'Box'
boxes = findobj(ax, 'Tag', 'Box');  % returns in reverse order (last to first)
boxes = flipud(boxes);              % now corresponds to positions 1..nCols

% map each column's group to a palette row
[~, gIdx] = ismember(string(go), string(cats));  % 1..K

for j = 1:nCols
    c = pal(gIdx(j),:);
    % Fill the box with a semi-transparent patch
    xd = get(boxes(j),'XData');
    yd = get(boxes(j),'YData');
    patch(ax, xd, yd, c, 'FaceAlpha', 0.35, 'EdgeColor', c, 'LineWidth', 1.2);
    % Hide the original box outline so the patch edge is the main edge
    set(boxes(j), 'LineStyle', 'none');
end

% Also color medians and whiskers (optional but usually nicer)
medLines = findobj(ax, 'Tag', 'Median');
medLines = flipud(medLines);
whiskers = findobj(ax, 'Tag', 'Whisker');
whiskers = flipud(whiskers);
caps = findobj(ax, 'Tag', 'Cap');
caps = flipud(caps);

for j = 1:nCols
    c = pal(gIdx(j),:);
    if j <= numel(medLines), set(medLines(j), 'Color', c, 'LineWidth', 1.5); end
end
% whiskers/caps: there are 2 per box (lower/upper)
for j = 1:nCols
    c = pal(gIdx(j),:);
    wIdx = (2*j-1):(2*j);
    if max(wIdx) <= numel(whiskers), set(whiskers(wIdx), 'Color', c, 'LineWidth', 1.0); end
    if max(wIdx) <= numel(caps),     set(caps(wIdx),     'Color', c, 'LineWidth', 1.0); end
end

% -------- legend (group colors)
% Build proxy patches for legend
proxy = gobjects(K,1);
for k = 1:K
    proxy(k) = patch(ax, nan, nan, pal(k,:), 'FaceAlpha', 0.35, 'EdgeColor', pal(k,:));
end
legend(ax, proxy, string(cats), 'Location', 'best');

hold(ax,'off');

% -------- output handles/info
h = struct();
h.boxplot = hBox;
h.order = ord;
h.sortedMedian = med(ord);
h.sortedGroup = go;
h.palette = pal;
h.groupCategories = cats;
end
