function [fig, stats] = create_wm_fig_cellfits(data_cell, data_cell2, data_cell3, Group_labels)
% CREATE_WM_FIG_CELLFITS Extract GM/WM/CSF values and plot mean cell-wise fits.
%
%   [fig, stats] = create_wm_fig_cellfits(data_cell, data_cell2, data_cell3, Group_labels)
%
% INPUTS
%   data_cell, data_cell2, data_cell3 - 1xS cell arrays of canlab objects
%   Group_labels                      - Sx1 labels used only for point colors
%
% OUTPUTS
%   fig   - figure handle
%   stats - struct with per-panel summaries
%
% NOTES
%   - Column 1 from extract_gray_white_csf is GM
%   - Column 2 is WM
%   - Column 3 is CSF
%   - Each cell/study is fit separately to obtain one slope/intercept
%   - Each panel shows one line using the mean cell-wise slope/intercept
%   - The shaded band reflects the 95% CI of the mean slope across cells

fun_wm = @(X) extract_gray_white_csf(X, 'masks', ...
    {'gray_matter_mask_sparse.img', 'canonical_white_matter.img', ...
     'canonical_ventricles.img'}, ...
    'eval', @(x1) median(x1, 1, 'omitnan'));

data_cell_wm1 = cellfun(fun_wm, data_cell, 'UniformOutput', false);
data_cell_wm2 = cellfun(fun_wm, data_cell2, 'UniformOutput', false);
data_cell_wm3 = cellfun(fun_wm, data_cell3, 'UniformOutput', false);

S = numel(data_cell_wm1);
if numel(data_cell_wm2) ~= S || numel(data_cell_wm3) ~= S || numel(Group_labels) ~= S
    error('All inputs must have the same number of cells/studies.');
end

groupCats = categorical(Group_labels);
groupLevels = categories(groupCats);
groupColors = lines(numel(groupLevels));
cellColors = zeros(S, 3);

for s = 1:S
    idx = find(strcmp(string(groupLevels), string(groupCats(s))), 1, 'first');
    cellColors(s, :) = groupColors(idx, :);
end

fig = figure('Position', [100 100 960 560]);

stats = struct();
stats.raw.csf = plot_panel(data_cell_wm1, 3, 1, cellColors, 'CSF', 'GM', 1, 'Raw: CSF vs GM');
stats.add.csf = plot_panel(data_cell_wm3, 3, 1, cellColors, 'CSF', 'GM', 2, 'Renormed (add only): CSF vs GM');
stats.addmul.csf = plot_panel(data_cell_wm2, 3, 1, cellColors, 'CSF', 'GM', 3, 'Renormed (add+mul): CSF vs GM');

stats.raw.wm = plot_panel(data_cell_wm1, 2, 1, cellColors, 'WM', 'GM', 4, 'Raw: WM vs GM');
stats.add.wm = plot_panel(data_cell_wm3, 2, 1, cellColors, 'WM', 'GM', 5, 'Renormed (add only): WM vs GM');
stats.addmul.wm = plot_panel(data_cell_wm2, 2, 1, cellColors, 'WM', 'GM', 6, 'Renormed (add+mul): WM vs GM');

end

function panelStats = plot_panel(dataCells, xcol, ycol, cellColors, xlabelStr, ylabelStr, subplotIdx, panelTitle)

subplot(2, 3, subplotIdx);
hold on

S = numel(dataCells);
slopes = nan(S, 1);
intercepts = nan(S, 1);
npts = nan(S, 1);
allx = [];
ally = [];
allColors = [];

for s = 1:S
    thisData = dataCells{s};
    if size(thisData, 2) < max(xcol, ycol)
        error('Each extracted cell must have at least %d columns.', max(xcol, ycol));
    end

    x = thisData(:, xcol);
    y = thisData(:, ycol);

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    if numel(x) < 2
        continue
    end

    p = polyfit(x, y, 1);
    slopes(s) = p(1);
    intercepts(s) = p(2);
    npts(s) = numel(x);

    allx = [allx; x]; %#ok<AGROW>
    ally = [ally; y]; %#ok<AGROW>
    allColors = [allColors; repmat(cellColors(s, :), numel(x), 1)]; %#ok<AGROW>
end

if isempty(allx)
    error('No valid data were available for panel: %s', panelTitle);
end

for i = 1:numel(allx)
    plot(allx(i), ally(i), '.', 'Color', allColors(i, :), 'MarkerSize', 10);
end

keep = isfinite(slopes) & isfinite(intercepts);
slopes = slopes(keep);
intercepts = intercepts(keep);
npts = npts(keep);

meanSlope = mean(slopes);
meanIntercept = mean(intercepts);

if numel(slopes) >= 2
    slopeSEM = std(slopes, 0) / sqrt(numel(slopes));
    slopeCI = meanSlope + tinv([0.025 0.975], numel(slopes) - 1) * slopeSEM;
else
    slopeCI = [NaN NaN];
end

xx = linspace(min(allx), max(allx), 200)';
yy = meanIntercept + meanSlope * xx;

if all(isfinite(slopeCI))
    yyLow = meanIntercept + slopeCI(1) * xx;
    yyHigh = meanIntercept + slopeCI(2) * xx;

    patch([xx; flipud(xx)], ...
          [yyLow; flipud(yyHigh)], ...
          [0.3 0.3 0.3], ...
          'FaceAlpha', 0.15, ...
          'EdgeColor', 'none');
end

plot(xx, yy, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2);

xlabel(xlabelStr);
ylabel(ylabelStr);
title(panelTitle);
box on
hold off

panelStats = struct();
panelStats.cellSlopes = slopes;
panelStats.cellIntercepts = intercepts;
panelStats.meanSlope = meanSlope;
panelStats.meanIntercept = meanIntercept;
panelStats.meanSlopeCI = slopeCI;
panelStats.nCells = numel(slopes);
panelStats.nPerCell = npts;

if numel(slopes) >= 2
    [h, p, ci, tstats] = ttest(slopes, 0);
    panelStats.ttest.h = h;
    panelStats.ttest.p = p;
    panelStats.ttest.ci = ci;
    panelStats.ttest.tstat = tstats.tstat;
    panelStats.ttest.df = tstats.df;
else
    panelStats.ttest = [];
end

end
