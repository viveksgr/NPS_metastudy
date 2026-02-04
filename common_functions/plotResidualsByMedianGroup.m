function [fig, stats] = plotResidualsByMedianGroup(x, y, varargin)
% plotResidualsByMedianGroup  Plot residuals split by median groups and t-test
%
% Usage:
%   [fig, stats] = plotResidualsByMedianGroup(x, y)
%   [fig, stats] = plotResidualsByMedianGroup(x, y, 'groupBy','both', 'nonparametric',true)
%
% Inputs:
%   x, y          - numeric vectors (same length). NaNs are removed pairwise.
%
% Optional name-value:
%   'groupBy'     - 'x' (default) | 'y' | 'both' | logical vector (same length as x)
%   'alpha'       - significance level (default 0.05)
%   'nonparametric' - true/false (default false). If true, uses ranksum (Wilcoxon)
%   'Jitter'      - scalar jitter amplitude for scatter (default 0.08)
%   'MarkerSize'  - scatter marker size (default 36)
%
% Outputs:
%   fig   - figure handle
%   stats - struct with fields:
%           .b         - regression betas [b0; b1]
%           .residuals - residual vector (same order after NaN-removal)
%           .groupIdx  - logical vector for group1 (true = group1)
%           .n1, .n2   - group sizes
%           .mean1,.mean2, .sem1, .sem2
%           .t, .p     - t-stat and p-value (or NaN if nonparametric requested)
%           .testName  - 'ttest2' or 'ranksum'
%
% Example:
%   rng(0); x = randn(200,1); y = 0.5*x + 0.2*randn(200,1);
%   plotResidualsByMedianGroup(x,y,'groupBy','both')

% ---------------- parse inputs --------------------------------
p = inputParser;
addParameter(p,'groupBy','x');
addParameter(p,'alpha',0.05,@(v)isnumeric(v)&&isscalar(v));
addParameter(p,'nonparametric',false,@islogical);
addParameter(p,'Jitter',0.08,@isnumeric);
addParameter(p,'MarkerSize',36,@isnumeric);
parse(p, varargin{:});
groupBy = p.Results.groupBy;
alpha = p.Results.alpha;
doRank = p.Results.nonparametric;
jitterAmp = p.Results.Jitter;
msz = p.Results.MarkerSize;

% ---------------- prepare data --------------------------------
if ~isvector(x) || ~isvector(y)
    error('x and y must be vectors.');
end
x = x(:); y = y(:);
valid = ~(isnan(x) | isnan(y));
x = x(valid); y = y(valid);
n = numel(x);
if n < 4
    error('Not enough observations after NaN removal (need >=4).');
end

% ---------------- fit linear regression y ~ x -------------------
X = [ones(n,1), x];
b = X \ y;                % OLS
yhat = X * b;
res = y - yhat;           % residuals (n x 1)

% ---------------- form grouping -------------------------------
if islogical(groupBy) && numel(groupBy) == n
    g1 = groupBy(:);
elseif ischar(groupBy) || isstring(groupBy)
    switch lower(string(groupBy))
        case 'x'
            med = median(x);
            g1 = x > med;                    % high-x group
        case 'y'
            med = median(y);
            g1 = y > med;                    % high-y group
        case 'both'
            gx = x > median(x);
            gy = y > median(y);
            g_highhigh = gx & gy;
            g_lowlow   = (~gx) & (~gy);
            % keep only high-high vs low-low; discard mixed rows
            keep = g_highhigh | g_lowlow;
            if sum(keep) < 4
                error('Not enough examples in high-high/low-low after median split.');
            end
            % restrict arrays to keep
            x = x(keep); y = y(keep); res = res(keep); n = numel(res);
            g1 = g_highhigh(keep);           % true = high-high
        otherwise
            error('Unknown ''groupBy'' option: %s', groupBy);
    end
else
    error('groupBy must be ''x'', ''y'', ''both'' or a logical vector of same length as inputs.');
end

% If groupBy was 'x' or 'y', we keep all rows and g1 already defined.
n1 = sum(g1);
n2 = sum(~g1);
if n1 < 2 || n2 < 2
    warning('One group has <2 members (n1=%d, n2=%d). t-test may be invalid.', n1, n2);
end

% ---------------- compute summary stats & test ------------------
r1 = res(g1);
r2 = res(~g1);
mean1 = mean(r1,'omitnan'); mean2 = mean(r2,'omitnan');
sem1  = std(r1,'omitnan')/sqrt(max(1,sum(~isnan(r1))));
sem2  = std(r2,'omitnan')/sqrt(max(1,sum(~isnan(r2))));

if doRank
    pval = ranksum(r1, r2);
    tstat = NaN;
    testName = 'ranksum';
else
    [h,pval,ci,stats_t] = ttest2(r1,r2,'Vartype','unequal','Alpha',alpha); % Welch
    tstat = stats_t.tstat;
    testName = 'ttest2 (Welch)';
end

% ---------------- make the bar plot -----------------------------
fig = figure('Color','w');
ax = axes('Parent',fig); hold(ax,'on');

% bar heights
barX = [1 2];
barHeights = [mean1 mean2];
barVals = bar(ax, barX, barHeights, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1.2);
% color: blue for group1, orange for group2
cols = [0.2 0.6 1; 1 .6 0.2];
for ii=1:2, barVals.CData(ii,:) = cols(ii,:); end

% error bars
errX = barX;
errorbar(ax, errX, barHeights, [sem1 sem2], 'k', 'LineStyle','none', 'LineWidth',1.2);

% jittered scatter of individual residuals
rng(0); % reproducible jitter
jx1 = barX(1) + (randn(numel(r1),1).*0.06) * jitterAmp;
jx2 = barX(2) + (randn(numel(r2),1).*0.06) * jitterAmp;
scatter(ax, jx1, r1, msz, repmat(cols(1,:), numel(r1),1), 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.8);
scatter(ax, jx2, r2, msz, repmat(cols(2,:), numel(r2),1), 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.8);

% aesthetics
plot(ax, xlim, [0 0], ':k'); % zero line
xlim([0.5 2.5]);
set(ax, 'XTick', barX, 'XTickLabel', {'Group 1','Group 2'}, 'TickLabelInterpreter','none');
ylabel('Residual (y - y_{pred})');
title(sprintf('Residuals by group (%s). %s: t=%.3g, p=%.3g', testName, testName, tstat, pval), 'Interpreter','none');

% annotate group Ns and means
text(1, barHeights(1)+1.6*sem1, sprintf('n=%d\nmean=%.3g', n1, mean1), 'HorizontalAlignment','center');
text(2, barHeights(2)+1.6*sem2, sprintf('n=%d\nmean=%.3g', n2, mean2), 'HorizontalAlignment','center');

% show p-value prominently
txt = sprintf('p = %.3g', pval);
if ~doRank && ~isnan(tstat)
    txt = sprintf('t = %.2f, p = %.3g', tstat, pval);
end
ylimVals = ylim;
text(mean(xlim), ylimVals(2) - 0.08*range(ylim), txt, 'HorizontalAlignment','center','FontWeight','bold');

box(ax,'on'); hold(ax,'off');

% ---------------- return stats struct --------------------------
stats = struct();
stats.b = b;
stats.residuals = res;
stats.groupIdx = g1;
stats.n1 = n1; stats.n2 = n2;
stats.mean1 = mean1; stats.mean2 = mean2;
stats.sem1 = sem1; stats.sem2 = sem2;
stats.t = tstat; stats.p = pval;
stats.testName = testName;

end
