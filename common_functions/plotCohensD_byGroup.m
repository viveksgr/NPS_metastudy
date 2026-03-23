function [fig, order] = plotCohensD_byGroup(outTbl, groupVec, varargin)
% PLOTCOHENSD_BYGROUP Plot Cohen's d (with SE) sorted by group means,
% with configurable horizontal spacing between groups.
%
% [fig, order] = plotCohensD_byGroup(outTbl, groupVec, 'GroupGap', 2, ...)
%
% See function body for options. New option:
%  'GroupGap' - integer >=0 number of blank slots inserted after each group (default 1)

% --- parse inputs
p = inputParser;
addParameter(p,'Colors',[], @(x) isnumeric(x) && (size(x,2)==3));
addParameter(p,'MarkerSize',36,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'XLabelRotation',45,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'ErrorBarWidth',0.4,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'FigureTitle','',@ischar);
addParameter(p,'Sorting',true,@islogical);
addParameter(p,'Axes',[], @(x) isempty(x) || isgraphics(x,'axes'));
addParameter(p,'ShowLegend',true,@islogical);
addParameter(p,'ShowXlabel',true,@islogical);
addParameter(p,'GroupGap',1,@(x)isnumeric(x) && isscalar(x) && (x>=0));
parse(p,varargin{:});

Cuser = p.Results.Colors;
ms = p.Results.MarkerSize;
rot = p.Results.XLabelRotation;
capw = p.Results.ErrorBarWidth;
figtitle = p.Results.FigureTitle;
sorting = p.Results.Sorting;
ax_in = p.Results.Axes;
showLegend = p.Results.ShowLegend;
showXlabel = p.Results.ShowXlabel;
groupGap = round(p.Results.GroupGap);

% --- sanity checks
nRows = height(outTbl);
if numel(groupVec) ~= nRows
    error('groupVec must be same length as number of rows in outTbl.');
end
if ~ismember('d', outTbl.Properties.VariableNames) || ~ismember('SE', outTbl.Properties.VariableNames)
    error('outTbl must contain columns named ''d'' and ''SE''.');
end

% extract names
if ismember('Name', outTbl.Properties.VariableNames)
    names = outTbl.Name;
else
    if ~isempty(outTbl.Properties.RowNames)
        names = outTbl.Properties.RowNames;
    else
        names = arrayfun(@(i) sprintf('Var%d', i), (1:nRows)', 'UniformOutput', false);
    end
end

% canonical group labels (for stable colour mapping)
grpCat = categorical(groupVec);
uniqueGroups = categories(grpCat);     % canonical label order
G = numel(uniqueGroups);

% prepare colours (one per unique group)
if isempty(Cuser)
    C = lines(G);  % default palette
else
    C = Cuser;
    if size(C,1) < G
        error('Provided Colors matrix has fewer rows (%d) than number of unique groups (%d).', size(C,1), G);
    end
    C = C(1:G, :);
end

% decide ordering (sorting on group mean + within-group d)
if sorting
    groupMeans = nan(G,1);
    for g = 1:G
        which = grpCat == uniqueGroups{g};
        groupMeans(g) = mean(outTbl.d(which), 'omitnan');
    end
    [~, grpOrderIdx] = sort(groupMeans, 'ascend');
    groups_sorted = uniqueGroups(grpOrderIdx);

    order = zeros(0,1);
    groupIdxLists = cell(G,1); % store original indices per group in display order
    for gi = 1:G
        gname = groups_sorted{gi};
        which = find(grpCat == gname);            % indices in original table
        [~, subord] = sort(outTbl.d(which), 'ascend', 'MissingPlacement','last');
        ord_this = which(subord);
        order = [order; ord_this(:)];
        groupIdxLists{gi} = ord_this(:);
    end
else
    order = (1:nRows).';
    groups_sorted = uniqueGroups;
    % preserve grouping order as they appear in groupVec unique order
    groupIdxLists = cell(G,1);
    for gi = 1:G
        groupIdxLists{gi} = find(grpCat == uniqueGroups{gi});
    end
end

% reorder table and names
Tsorted = outTbl(order, :);
names_sorted = names(order);
grp_sorted = grpCat(order);

% compute x positions with gaps between groups
xpos = zeros(numel(order),1);
current = 1;
groupBoundaries = zeros(G,2); % start,end xpos in the final axis for separators
idx = 1;
for gi = 1:G
    members = groupIdxLists{gi};
    nmem = numel(members);
    if nmem == 0
        groupBoundaries(gi,:) = [NaN NaN];
        continue;
    end
    % assign consecutive positions for this group
    xpos(idx:idx+nmem-1) = current:(current + nmem - 1);
    startpos = current;
    endpos = current + nmem - 1;
    groupBoundaries(gi,:) = [startpos, endpos];
    % advance current by group size + gap
    current = endpos + 1 + groupGap;
    idx = idx + nmem;
end

% map sorted rows to canonical group colour index (stable mapping)
colorIdx = zeros(numel(order),1);
for i = 1:numel(order)
    colorIdx(i) = find(strcmp(string(uniqueGroups), string(grp_sorted(i))));
end

% --- plotting
ax = ax_in;
if isempty(ax)
    fig = figure('Color','w','Position',[100 200 1100 420]);
    ax = axes('Parent',fig);
else
    fig = ancestor(ax,'figure');
end
hold(ax,'on');

x = xpos;
y = Tsorted.d;
se = Tsorted.SE;

% Plot errorbars (custom caps) at custom x positions
for i = 1:numel(x)
    xi = x(i);
    yi = y(i);
    s = se(i)*1.96;
    % vertical line
    line(ax, [xi xi], [yi - s, yi + s], 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
    % horizontal caps
    cap = capw * 0.5; % half-width (in axis units)
    line(ax, [xi-cap xi+cap], [yi - s yi - s], 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
    line(ax, [xi-cap xi+cap], [yi + s yi + s], 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
end

% scatter points colored by canonical-group colors (stable mapping)
for g = 1:G
    wh = find(colorIdx == g);
    if ~isempty(wh)
        scatter(ax, x(wh), y(wh), ms, repmat(C(g,:), numel(wh), 1), 'filled', 'MarkerEdgeColor', 'k');
    end
end

% subtle separators between groups (draw in background)
yl = ylim(gca); % use current y limits to draw lines, will update after setting y-limits more precisely
for gi = 1:(G-1)
    b = groupBoundaries(gi,2);
    if ~isnan(b)
        sepx = b + groupGap/2;
        % draw faint vertical line if gap>0
        % if groupGap > 0
        %     line(ax, [sepx sepx], [yl(1) yl(2)], 'Color', [0.85 0.85 0.85], 'LineStyle','-', 'LineWidth', 0.8, 'HandleVisibility','off');
        % end
    end
end

% Plot zero line
ymin = min(y - se);
ymax = max(y + se);
if ymin==ymax
    ymin = ymin - 1; ymax = ymax + 1;
end
ylimPadding = 0.08 * (ymax - ymin + eps);
ylim(ax, [ymin - ylimPadding, ymax + ylimPadding]);
plot(ax, get(ax,'XLim'), [0 0], ':k', 'LineWidth', 1);

% X ticks/labels at the point locations
if showXlabel
    set(ax, 'XTick', x, 'XTickLabel', names_sorted, 'TickLabelInterpreter', 'none');
    xtickangle(rot);
else
    set(ax, 'XTick', [], 'XTickLabel', [], 'TickLabelInterpreter', 'none');
end

% labels and title
ylabel(ax, 'Cohen''s d');
if ~isempty(figtitle)
    title(ax, figtitle, 'Interpreter', 'none');
end

% Legend: use canonical group order (so legend colours are stable across calls)
legendHandles = gobjects(G,1);
legendLabels = cell(G,1);
for g = 1:G
    legendHandles(g) = scatter(ax, NaN, NaN, ms, C(g,:), 'filled', 'MarkerEdgeColor', 'k');
    legendLabels{g} = sprintf('%s', string(uniqueGroups{g}));
end
if showLegend
    legend(ax, legendHandles, legendLabels, 'Location', 'bestoutside');
end

box(ax,'on');
grid(ax,'off');
hold(ax,'off');

% return order as column vector
order = order(:);

end
