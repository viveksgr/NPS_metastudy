function plotTest0ByGroup(R, cl_mat, ttl, varargin)
% plotTest0ByGroup  Bar plot of per-group Test-0 decoding accuracy, with the
%                   permutation floor and pooled ceiling drawn as lines on
%                   each bar and significance shown as asterisks.
%
%   plotTest0ByGroup(R, cl_mat, 'Test 0 - random CV')
%   plotTest0ByGroup(R, cl_mat, ttl, 'sortBars', true)
%
%   R       - table from run_test0_bygroup (Group, Acc, Floor, Ceiling, p, N)
%   cl_mat  - Nx3 RGB, one row per group in the ORIGINAL (unsorted) R order.
%   ttl     - title string.
%
%   Name-value options:
%     'sortBars' - false (default) keeps the native group order so that the
%                  CV and LOSO figures are directly comparable bar-for-bar.
%                  true sorts ascending by accuracy (plotFixedEffects style).
%     'yLabel'   - y-axis label.
%
%   Styling follows common_functions/plotFixedEffects.m.

    p = inputParser;
    addParameter(p, 'sortBars', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'yLabel', 'Decoding accuracy (per-group recall)', @(x) ischar(x)||isstring(x));
    parse(p, varargin{:});
    sortBars = p.Results.sortBars;

    acc   = R.Acc;
    flo   = R.Floor;
    ceil_ = R.Ceiling;
    pval  = R.p;
    names = string(R.Group);
    n     = numel(acc);

    if size(cl_mat,1) < n
        error('cl_mat must have at least %d rows (one per group).', n);
    end

    if sortBars
        [~, idx] = sort(acc, 'ascend'); idx = idx(:);
        acc = acc(idx); flo = flo(idx); ceil_ = ceil_(idx);
        pval = pval(idx); names = names(idx); cl_mat = cl_mat(idx,:);
    end

    figure('Position',[0.5 0.5 460 300]); hold on;

    for i = 1:n
        bar(i, acc(i), 'FaceColor', cl_mat(i,:), 'EdgeColor', 'k');
    end

    % floor (permutation null) and ceiling (pooled, matched N) as lines on bars
    hw = 0.42;
    hF = gobjects(1); hC = gobjects(1);
    for i = 1:n
        if ~isnan(flo(i))
            h = plot([i-hw i+hw], [flo(i) flo(i)], 'k--', 'LineWidth', 1.3);
            if i==1, hF = h; end
        end
        if ~isnan(ceil_(i))
            h = plot([i-hw i+hw], [ceil_(i) ceil_(i)], '-', ...
                     'Color', [0.75 0 0], 'LineWidth', 1.6);
            if i==1, hC = h; end
        end
    end

    allv = [acc; flo; ceil_];
    yr = max(allv,[],'omitnan') - min([0; min(allv,[],'omitnan')]);
    offset = 0.03 * yr;
    for i = 1:n
        stars = pStars(pval(i));
        if isempty(stars), continue; end
        top = max([acc(i), ceil_(i)], [], 'omitnan');
        text(i, top + offset, stars, 'Color', 'r', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontWeight', 'bold');
    end

    set(gca, 'XTick', 1:n, 'XTickLabel', cellstr(names), ...
             'XTickLabelRotation', 45, ...
             'TickLabelInterpreter', 'none');
    xlim([0.4 n+0.6]);
    ylim([0, max(allv,[],'omitnan') + 0.12*yr]);
    ylabel(p.Results.yLabel);
    title(ttl, 'Interpreter','none');
    if isgraphics(hF) && isgraphics(hC)
        legend([hF hC], {'permutation floor','pooled ceiling (matched N)'}, ...
               'Location','northoutside','Orientation','horizontal','Box','off');
    end
    box off;
    hold off;
end

function s = pStars(p)
    if     isnan(p),  s = '';
    elseif p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = '';
    end
end
