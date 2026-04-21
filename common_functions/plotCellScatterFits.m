function stats = plotCellScatterFits(C1, C2, scatterColor)
% plotCellScatterMeanSlope
% Plot all cell-pair points in one color, compute one slope from the
% individual cell slopes, and show that overall slope with a shaded band.
%
% Usage:
%   stats = plotCellScatterMeanSlope(C1, C2)
%   stats = plotCellScatterMeanSlope(C1, C2, [0 0.45 0.74])
%
% Inputs:
%   C1, C2       same-sized cell arrays of numeric vectors
%   scatterColor RGB triplet for all points and the overall fit line
%
% Output:
%   stats        struct with per-cell slopes and overall slope stats
%
% Overall slope:
%   This function first fits each cell separately:
%       y_i = m_i*x + b_i
%   Then computes one overall slope as:
%       mean(m_i)
%   and tests whether the cell slopes differ from 0.

    if nargin < 3
        scatterColor = [0 0.45 0.74];
    end

    if ~iscell(C1) || ~iscell(C2)
        error('C1 and C2 must both be cell arrays.');
    end

    if ~isequal(size(C1), size(C2))
        error('C1 and C2 must have the same size.');
    end

    n = numel(C1);
    slopes = nan(n,1);
    intercepts = nan(n,1);
    allx = [];
    ally = [];

    for i = 1:n
        x = C1{i};
        y = C2{i};

        if ~isvector(x) || ~isvector(y)
            error('Each entry of C1 and C2 must be a vector. Problem at cell %d.', i);
        end

        x = x(:);
        y = y(:);

        if numel(x) ~= numel(y)
            error('C1{%d} and C2{%d} must have equal length.', i, i);
        end

        valid = isfinite(x) & isfinite(y);
        x = x(valid);
        y = y(valid);

        if numel(x) < 2
            warning('Cell %d has fewer than 2 valid points. Skipping.', i);
            continue
        end

        p = polyfit(x, y, 1);
        slopes(i) = p(1);
        intercepts(i) = p(2);

        allx = [allx; x];
        ally = [ally; y];
    end

    keep = isfinite(slopes);
    slopes = slopes(keep);
    intercepts = intercepts(keep);

    if isempty(slopes)
        error('No valid cell fits could be computed.');
    end

    % One overall slope from the individual cell slopes
    meanSlope = mean(slopes);

    % Use the mean intercept too, so the displayed line is based on the
    % cell-wise fits rather than a pooled fit to all raw points.
    meanIntercept = mean(intercepts);

    % Stats on slopes across cells
    stats = struct();
    stats.cellSlopes = slopes;
    stats.cellIntercepts = intercepts;
    stats.meanSlope = meanSlope;
    stats.meanIntercept = meanIntercept;
    stats.nCells = numel(slopes);

    if numel(slopes) >= 2
        [h,p,ci,ttestStats] = ttest(slopes, 0);
        stats.ttest.h = h;
        stats.ttest.p = p;
        stats.ttest.ci = ci;
        stats.ttest.tstat = ttestStats.tstat;
        stats.ttest.df = ttestStats.df;

        stats.signrank.p = signrank(slopes, 0);

        slopeSEM = std(slopes) / sqrt(numel(slopes));
        slopeCI = tinv([0.025 0.975], numel(slopes)-1) * slopeSEM + meanSlope;
        stats.meanSlopeCI = slopeCI;
    else
        stats.ttest = [];
        stats.signrank = [];
        stats.meanSlopeCI = [NaN NaN];
    end

    % Plot
    figure('Position',[0.5 0.5 320 240]);
    hold on

    plot(allx, ally, '.', 'Color', scatterColor, 'MarkerSize', 10);

    xx = linspace(min(allx), max(allx), 200)';
    yy = meanIntercept + meanSlope * xx;

    if all(isfinite(stats.meanSlopeCI))
        yyLow = meanIntercept + stats.meanSlopeCI(1) * xx;
        yyHigh = meanIntercept + stats.meanSlopeCI(2) * xx;

        patch([xx; flipud(xx)], ...
              [yyLow; flipud(yyHigh)], ...
              scatterColor, ...
              'FaceAlpha', 0.15, ...
              'EdgeColor', 'none');
    end

    plot(xx, yy, '-', 'Color', scatterColor, 'LineWidth', 2);

    xlabel('C1');
    ylabel('C2');
    title('All Points with One Overall Slope from Cell-wise Slopes');
    box on
    hold off

    % Command window summary
    fprintf('\nOverall slope from cell-wise slopes\n');
    fprintf('Number of cells used: %d\n', stats.nCells);
    fprintf('Mean slope: %.6f\n', stats.meanSlope);

    if ~isempty(stats.ttest)
        fprintf('t-test vs 0: t(%d) = %.4f, p = %.4g\n', ...
            stats.ttest.df, stats.ttest.tstat, stats.ttest.p);
        fprintf('95%% CI for mean slope: [%.6f, %.6f]\n', ...
            stats.ttest.ci(1), stats.ttest.ci(2));
        fprintf('Sign-rank p = %.4g\n', stats.signrank.p);
    end
end
