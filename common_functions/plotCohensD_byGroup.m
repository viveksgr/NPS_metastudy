function [fig, order] = plotCohensD_byGroup(outTbl, groupVec, varargin)
% PLOTCOHENSD_BYGROUP Plot Cohen's d (with SE) sorted by group means.
% [fig, order] = plotCohensD_byGroup(outTbl, groupVec, 'Sorting', true, 'Colors', C, ...)
%
% Key changes vs your previous version:
% - Colors are tied to group *labels* (unique groups) so colors remain fixed across runs.
% - Returns 'order' used to reorder rows. When Sorting==false, order = (1:nRows)'.

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

parse(p,varargin{:});
Cuser = p.Results.Colors;
ms = p.Results.MarkerSize;
rot = p.Results.XLabelRotation;
capw = p.Results.ErrorBarWidth;
figtitle = p.Results.FigureTitle;
sorting = p.Results.Sorting;

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

% canonical group labels (used for stable color mapping)
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
    % if user provided more rows than G, we only use first G rows
    C = C(1:G, :);
end

% decide ordering (sorting on group mean + within-group d)
if sorting
    % compute group means (over outTbl.d) for canonical group labels
    groupMeans = nan(G,1);
    for g = 1:G
        which = grpCat == uniqueGroups{g};
        groupMeans(g) = mean(outTbl.d(which), 'omitnan');
    end
    % sort groups by their mean (ascending) to produce display order
    [~, grpOrderIdx] = sort(groupMeans, 'ascend');
    groups_sorted = uniqueGroups(grpOrderIdx);

    % now build order by concatenating within-group sorted-by-d indices
    order = zeros(0,1);
    for gi = 1:G
        gname = groups_sorted{gi};
        which = find(grpCat == gname);            % indices in original table
        % sort these indices by their d ascending
        [~, subord] = sort(outTbl.d(which), 'ascend', 'MissingPlacement','last');
        ord_this = which(subord);
        order = [order; ord_this(:)];
    end
else
    % keep original order
    order = (1:nRows).';
    groups_sorted = uniqueGroups; % not used for ordering but keep defined
end

% reorder table and names according to order
Tsorted = outTbl(order, :);
names_sorted = names(order);
grp_sorted = grpCat(order);

% map each sorted row to the colour index based on canonical label order
colorIdx = zeros(numel(order),1);
for i = 1:numel(order)
    % find which canonical group this row belongs to (stable mapping)
    colorIdx(i) = find(strcmp(string(uniqueGroups), string(grp_sorted(i))));
end

% --- plotting
ax = p.Results.Axes;
if isempty(ax)
    fig = figure('Color','w','Position',[100 200 1100 420]);
    ax = axes('Parent',fig);
    createdFigure = true;
else
    % do not create a figure; use supplied axes
    fig = ancestor(ax,'figure');
    createdFigure = false;
end
hold(ax,'on');


x = 1:height(Tsorted);
y = Tsorted.d;
se = Tsorted.SE;

% Plot errorbars (custom caps)
for i = 1:numel(x)
    xi = x(i);
    yi = y(i);
    s = se(i);
    % vertical line
    line([xi xi], [yi - s, yi + s], 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
    % horizontal caps
    cap = capw * 0.5; % half-width
    line([xi-cap xi+cap], [yi - s yi - s], 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
    line([xi-cap xi+cap], [yi + s yi + s], 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
end

% scatter points colored by canonical-group colors (stable mapping)
for g = 1:G
    wh = find(colorIdx == g);
    if ~isempty(wh)
        scatter(ax, x(wh), y(wh), ms, repmat(C(g,:), numel(wh), 1), 'filled', 'MarkerEdgeColor', 'k');
    end
end

% Plot zero line
ymin = min(y - se);
ymax = max(y + se);
ylimPadding = 0.08 * (ymax - ymin + eps);
ylim([ymin - ylimPadding, ymax + ylimPadding]);
plot(xlim, [0 0], ':k', 'LineWidth', 1);

% X ticks/labels
if  p.Results.ShowXlabel
    set(ax, 'XTick', x, 'XTickLabel', names_sorted, 'TickLabelInterpreter', 'none');
    xtickangle(rot);
else
    set(ax, 'XTick', [], 'XTickLabel', [], 'TickLabelInterpreter', 'none');
end

% labels and title
ylabel('Cohen''s d');
if ~isempty(figtitle)
    title(figtitle, 'Interpreter', 'none');
end

% Legend: use canonical group order (so legend colors are stable across calls)
legendHandles = gobjects(G,1);
legendLabels = cell(G,1);
for g = 1:G
    legendHandles(g) = scatter(NaN, NaN, ms, C(g,:), 'filled', 'MarkerEdgeColor', 'k');
    legendLabels{g} = sprintf('%s', string(uniqueGroups{g}));
end
if p.Results.ShowLegend
    legend(legendHandles, legendLabels, 'Location', 'bestoutside');
end


box(ax,'on');
grid(ax,'off');
hold(ax,'off');

% return order as column vector
order = order(:);

end
