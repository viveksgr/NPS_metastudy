function [fig, ax, stats] = plotCohensD_bivariate_byGroup(outTbl1, outTbl2, groupVec, varargin)
% PLOTCOHENSD_BIVARIATE_BYGROUP  Scatter Cohen's d from two tables against
% each other, with 95% CI error bars on both axes, colored by group.
%
%   [fig, ax] = plotCohensD_bivariate_byGroup(outTbl1, outTbl2, groupVec, ...)
%
% INPUTS
%   outTbl1  - table with columns: Name, d, SE (n optional). Plotted on x.
%   outTbl2  - table with columns: Name, d, SE (n optional). Plotted on y.
%              Must have same number of rows as outTbl1; rows are assumed
%              to correspond 1:1. A warning is issued if Name columns
%              disagree.
%   groupVec - group label per row (numeric, string, cell, or categorical)
%
% NAME/VALUE OPTIONS
%   'Colors'        - G x 3 color matrix (default: lines(G))
%   'MarkerSize'    - marker area (default 50)
%   'XLabel'        - x-axis label (default 'd_1')
%   'YLabel'        - y-axis label (default 'd_2')
%   'FigureTitle'   - figure title (default '')
%   'Axes'          - target axes handle (default: new figure)
%   'ShowLegend'    - true/false (default true)
%   'ShowLabels'    - annotate each point with its Name (default false)
%   'UnityLine'     - draw y=x reference line (default false)
%   'RegressionLine' - draw OLS line of y on x with slope test
%                      (default true). Replaces the unity line visually.
%   'ShowStatsBox'  - show slope, p, R^2 annotation on plot (default true)
%   'ZeroLines'     - draw x=0 and y=0 reference lines (default true)
%   'AxisEqual'     - force equal x/y scaling (default true)
%   'GroupRegions'  - draw faded colored region per group (default false)
%   'GroupRegionStyle' - 'hull' (convex hull of points) or 'ellipse'
%                        (covariance-based ~95% coverage). Default 'hull'.
%   'GroupRegionAlpha' - face transparency for group regions (default 0.15)
%   'XLim'          - [xmin xmax], overrides auto-computed x-limits
%   'YLim'          - [ymin ymax], overrides auto-computed y-limits
%
% OUTPUTS
%   fig   - figure handle
%   ax    - axes handle
%   stats - struct with fields: slope, intercept, slope_SE, slope_p,
%           slope_CI95, R2, n (based on OLS regression of y on x)

% --- parse inputs
p = inputParser;
addParameter(p,'Colors',[], @(x) isnumeric(x) && (size(x,2)==3));
addParameter(p,'MarkerSize',50,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'XLabel','d_1',@(x)ischar(x)||isstring(x));
addParameter(p,'YLabel','d_2',@(x)ischar(x)||isstring(x));
addParameter(p,'FigureTitle','',@(x)ischar(x)||isstring(x));
addParameter(p,'Axes',[], @(x) isempty(x) || isgraphics(x,'axes'));
addParameter(p,'ShowLegend',true,@islogical);
addParameter(p,'ShowLabels',false,@islogical);
addParameter(p,'UnityLine',false,@islogical);
addParameter(p,'RegressionLine',true,@islogical);
addParameter(p,'ShowStatsBox',true,@islogical);
addParameter(p,'ZeroLines',true,@islogical);
addParameter(p,'AxisEqual',true,@islogical);
addParameter(p,'GroupRegions',false,@islogical);
addParameter(p,'GroupRegionStyle','hull', ...
    @(s) any(strcmpi(s,{'hull','ellipse'})));
addParameter(p,'GroupRegionAlpha',0.15, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
addParameter(p,'XLim',[], ...
    @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
addParameter(p,'YLim',[], ...
    @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
parse(p,varargin{:});

Cuser       = p.Results.Colors;
ms          = p.Results.MarkerSize;
xlab        = char(p.Results.XLabel);
ylab        = char(p.Results.YLabel);
figtitle    = char(p.Results.FigureTitle);
ax_in       = p.Results.Axes;
showLegend  = p.Results.ShowLegend;
showLabels  = p.Results.ShowLabels;
doUnity     = p.Results.UnityLine;
doReg       = p.Results.RegressionLine;
doStatsBox  = p.Results.ShowStatsBox;
doZero      = p.Results.ZeroLines;
axEqual     = p.Results.AxisEqual;
doRegions   = p.Results.GroupRegions;
regStyle    = lower(p.Results.GroupRegionStyle);
regAlpha    = p.Results.GroupRegionAlpha;
xlimUser    = p.Results.XLim;
ylimUser    = p.Results.YLim;

% --- sanity checks
nRows = height(outTbl1);
if height(outTbl2) ~= nRows
    error('outTbl1 and outTbl2 must have the same number of rows.');
end
if numel(groupVec) ~= nRows
    error('groupVec must have same length as rows in the tables.');
end
required = {'d','SE'};
miss1 = setdiff(required, outTbl1.Properties.VariableNames);
miss2 = setdiff(required, outTbl2.Properties.VariableNames);
if ~isempty(miss1)
    error('outTbl1 missing required column(s): %s', strjoin(miss1, ', '));
end
if ~isempty(miss2)
    error('outTbl2 missing required column(s): %s', strjoin(miss2, ', '));
end

% extract names from outTbl1 (preferred); warn if outTbl2 disagrees
if ismember('Name', outTbl1.Properties.VariableNames)
    names = outTbl1.Name;
elseif ~isempty(outTbl1.Properties.RowNames)
    names = outTbl1.Properties.RowNames;
else
    names = arrayfun(@(i) sprintf('Var%d', i), (1:nRows)', 'UniformOutput', false);
end
if ismember('Name', outTbl2.Properties.VariableNames)
    names2 = outTbl2.Name;
    if ~isequal(string(names(:)), string(names2(:)))
        warning('plotCohensD_bivariate_byGroup:NameMismatch', ...
            'Name columns differ between outTbl1 and outTbl2. Rows are matched by position.');
    end
end

% canonical groups
grpCat = categorical(groupVec);
uniqueGroups = categories(grpCat);
G = numel(uniqueGroups);

% colors
if isempty(Cuser)
    C = lines(G);
else
    if size(Cuser,1) < G
        error('Colors has fewer rows (%d) than groups (%d).', size(Cuser,1), G);
    end
    C = Cuser(1:G, :);
end

% pull data (x from table1, y from table2)
x  = outTbl1.d;
y  = outTbl2.d;
xe = outTbl1.SE;
ye = outTbl2.SE;

% --- OLS regression of y on x (for line + stats panel)
stats = struct('slope',NaN,'intercept',NaN,'slope_SE',NaN,'slope_p',NaN, ...
               'slope_CI95',[NaN NaN],'R2',NaN,'n',0);
valid = isfinite(x) & isfinite(y);
nValid = sum(valid);
if nValid >= 3
    xv = x(valid);
    yv = y(valid);
    Xd = [ones(nValid,1), xv];
    b  = Xd \ yv;                              % [intercept; slope]
    yhat = Xd * b;
    resid = yv - yhat;
    dfR = nValid - 2;
    sigma2 = sum(resid.^2) / dfR;
    covB   = sigma2 * inv(Xd' * Xd);           %#ok<MINV>
    se_b1  = sqrt(covB(2,2));
    t_b1   = b(2) / se_b1;
    p_b1   = 2 * tcdf(-abs(t_b1), dfR);
    tcrit  = tinv(0.975, dfR);
    ci_b1  = [b(2) - tcrit*se_b1, b(2) + tcrit*se_b1];
    SStot  = sum((yv - mean(yv)).^2);
    SSres  = sum(resid.^2);
    R2     = 1 - SSres/SStot;

    stats.slope      = b(2);
    stats.intercept  = b(1);
    stats.slope_SE   = se_b1;
    stats.slope_p    = p_b1;
    stats.slope_CI95 = ci_b1;
    stats.R2         = R2;
    stats.n          = nValid;
end

% --- axes setup
if isempty(ax_in)
    fig = figure('Color','w','Position',[100 200 620 560]);
    ax  = axes('Parent',fig);
else
    ax  = ax_in;
    fig = ancestor(ax,'figure');
end
hold(ax,'on');

% --- optional group regions (drawn first, underneath everything)
if doRegions
    for g = 1:G
        wh = grpCat == uniqueGroups{g};
        gx = x(wh);
        gy = y(wh);
        valid = isfinite(gx) & isfinite(gy);
        gx = gx(valid);
        gy = gy(valid);
        np = numel(gx);
        if np < 1, continue; end

        col = C(g,:);

        switch regStyle
            case 'hull'
                if np >= 3
                    k = convhull(gx, gy);
                    patch(ax, gx(k), gy(k), col, ...
                        'FaceAlpha', regAlpha, ...
                        'EdgeColor', col, 'EdgeAlpha', regAlpha*2, ...
                        'LineWidth', 0.75, ...
                        'HandleVisibility','off');
                else
                    % fewer than 3 points: draw a small disc centered on mean
                    cx = mean(gx); cy = mean(gy);
                    r  = max([std(gx), std(gy), 0.05], [], 'omitnan');
                    tt = linspace(0, 2*pi, 60);
                    patch(ax, cx + r*cos(tt), cy + r*sin(tt), col, ...
                        'FaceAlpha', regAlpha, ...
                        'EdgeColor', col, 'EdgeAlpha', regAlpha*2, ...
                        'LineWidth', 0.75, ...
                        'HandleVisibility','off');
                end

            case 'ellipse'
                % covariance ellipse at ~95% for bivariate normal
                % (scale factor sqrt(5.991) ~ chi2inv(0.95, 2))
                if np >= 2
                    mu = [mean(gx), mean(gy)];
                    Sg = cov(gx, gy);
                    % handle degenerate covariance
                    if any(~isfinite(Sg(:))) || rank(Sg) < 2
                        r = max([std(gx), std(gy), 0.05], [], 'omitnan');
                        tt = linspace(0, 2*pi, 60);
                        ex = mu(1) + r*cos(tt);
                        ey = mu(2) + r*sin(tt);
                    else
                        [V, D] = eig(Sg);
                        tt = linspace(0, 2*pi, 120);
                        circ = [cos(tt); sin(tt)];
                        ell  = V * sqrt(D) * circ * sqrt(5.991);
                        ex = mu(1) + ell(1,:);
                        ey = mu(2) + ell(2,:);
                    end
                    patch(ax, ex, ey, col, ...
                        'FaceAlpha', regAlpha, ...
                        'EdgeColor', col, 'EdgeAlpha', regAlpha*2, ...
                        'LineWidth', 0.75, ...
                        'HandleVisibility','off');
                else
                    % single point: small disc fallback
                    cx = gx(1); cy = gy(1);
                    r  = 0.05;
                    tt = linspace(0, 2*pi, 60);
                    patch(ax, cx + r*cos(tt), cy + r*sin(tt), col, ...
                        'FaceAlpha', regAlpha, ...
                        'EdgeColor', col, 'EdgeAlpha', regAlpha*2, ...
                        'LineWidth', 0.75, ...
                        'HandleVisibility','off');
                end
        end
    end
end

% --- error bars (drawn first, underneath points)
errAlpha = 0.35;
for i = 1:nRows
    % find group color for this row
    g = find(strcmp(string(uniqueGroups), string(grpCat(i))), 1);
    col = [C(g,:), errAlpha];
    % horizontal bar (x error)
    line(ax, [x(i)-xe(i), x(i)+xe(i)], [y(i) y(i)], ...
        'Color', col, 'LineWidth', 1);
    % vertical bar (y error)
    line(ax, [x(i) x(i)], [y(i)-ye(i), y(i)+ye(i)], ...
        'Color', col, 'LineWidth', 1);
end

% --- scatter points per group (for a clean legend)
hGroup = gobjects(G,1);
for g = 1:G
    wh = grpCat == uniqueGroups{g};
    hGroup(g) = scatter(ax, x(wh), y(wh), ms, C(g,:), 'filled', ...
        'MarkerEdgeColor','k');
end

% --- reference lines
xl = [min(x - xe), max(x + xe)];
yl = [min(y - ye), max(y + ye)];
pad = 0.08 * max(diff(xl), diff(yl));
xl = xl + [-pad, pad];
yl = yl + [-pad, pad];

if axEqual
    lo = min(xl(1), yl(1));
    hi = max(xl(2), yl(2));
    xl = [lo hi];
    yl = [lo hi];
end

% user-provided limits override auto-computed ones
if ~isempty(xlimUser), xl = xlimUser(:)'; end
if ~isempty(ylimUser), yl = ylimUser(:)'; end

xlim(ax, xl);
ylim(ax, yl);

if doZero
    plot(ax, xl, [0 0], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility','off');
    plot(ax, [0 0], yl, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility','off');
end

if doUnity
    lo = max(xl(1), yl(1));
    hi = min(xl(2), yl(2));
    plot(ax, [lo hi], [lo hi], ':', 'Color', [0.6 0.6 0.6], ...
        'LineWidth', 1, 'HandleVisibility','off');
end

if doReg && isfinite(stats.slope)
    xx = linspace(xl(1), xl(2), 100);
    yy = stats.intercept + stats.slope * xx;
    % solid black if significant, dashed gray if not
    if stats.slope_p < 0.05
        plot(ax, xx, yy, '-', 'Color', [0.1 0.1 0.1], ...
            'LineWidth', 1.5, 'HandleVisibility','off');
    else
        plot(ax, xx, yy, '--', 'Color', [0.4 0.4 0.4], ...
            'LineWidth', 1.2, 'HandleVisibility','off');
    end
end

% --- stats annotation box
if doStatsBox && isfinite(stats.slope)
    if stats.slope_p < 1e-3
        pstr = 'p < 0.001';
    else
        pstr = sprintf('p = %.3f', stats.slope_p);
    end
    if stats.slope_p < 0.05, sig_mark = ' *'; else, sig_mark = ''; end
    txt = sprintf('slope = %.3f \\pm %.3f%s\n%s\nR^2 = %.3f   n = %d', ...
        stats.slope, stats.slope_SE, sig_mark, pstr, stats.R2, stats.n);
    text(ax, 0.03, 0.97, txt, 'Units','normalized', ...
        'VerticalAlignment','top', 'HorizontalAlignment','left', ...
        'FontSize', 9, 'BackgroundColor', [1 1 1 0.75], ...
        'EdgeColor', [0.7 0.7 0.7], 'Margin', 4);
end

if axEqual
    % axis(ax, 'equal');
    xlim(ax, xl);
    ylim(ax, yl);
end

% --- optional point labels
if showLabels
    for i = 1:nRows
        text(ax, x(i), y(i), ['  ' char(string(names{i}))], ...
            'FontSize', 8, 'Interpreter', 'none');
    end
end

% --- cosmetics
xlabel(ax, xlab);
ylabel(ax, ylab);
if ~isempty(figtitle)
    title(ax, figtitle, 'Interpreter', 'none');
end
box(ax,'on');
grid(ax,'off');

if showLegend
    legend(ax, hGroup, cellstr(uniqueGroups), 'Location', 'bestoutside');
end

hold(ax,'off');

end