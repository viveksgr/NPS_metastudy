function plotFixedEffects(stats, cl_mat2, varargin)
% plotFixedEffects  Bar plot of fixed-effect betas with SE error bars,
%                   sorted by descending beta.
%
%   plotFixedEffects(stats, cl_mat2)
%   plotFixedEffects(stats, cl_mat2, 'specialIntercept', false)
%
%   stats   - struct/dataset returned by fixedEffects (must have .Estimate,
%             .SE, and .Name fields/columns).
%   cl_mat2 - Nx3 matrix of RGB colors, one row per coefficient (in the
%             original, unsorted order from stats).
%
%   Name-value options:
%     'specialIntercept' - if true (default), keep the first coefficient in
%                          first position and color it gray. If false, sort
%                          all coefficients together. Colors are always
%                          reordered with their original coefficient rows.

    p = inputParser;
    addParameter(p, 'specialIntercept', true, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    specialIntercept = p.Results.specialIntercept;

    betas = stats.Estimate;
    lo    = stats.Lower;
    up    = stats.Upper;
    pval  = stats.pValue;
    names = stats.Name;
    n     = numel(betas);

    if size(cl_mat2,1) < n
        error('cl_mat2 must have at least %d rows (one per coefficient).', n);
    end

    originalIdx = (1:n)';
    if specialIntercept
        [~, sortIdx] = sort(betas(2:end), 'ascend');
        idx = [1; sortIdx(:) + 1];
    else
        [~, idx] = sort(betas, 'ascend');
        idx = idx(:);
    end

    betas   = betas(idx);
    lo      = lo(idx);
    up      = up(idx);
    pval    = pval(idx);
    names   = names(idx);
    cl_mat2 = cl_mat2(idx,:);
    originalIdx = originalIdx(idx);

    if specialIntercept
        cl_mat2(originalIdx == 1,:) = [0.5 0.5 0.5];
    end

    negLen = betas - lo;
    posLen = up - betas;

    figure('Position',[0.5 0.5 320 240]); hold on;
    for i = 1:n
        bar(i, betas(i), 'FaceColor', cl_mat2(i,:), 'EdgeColor', 'k');
    end
    errorbar(1:n, betas, negLen, posLen, 'k', 'linestyle', 'none', 'LineWidth', 1.2, 'CapSize', 8);

    yr     = max(up) - min(lo);
    offset = 0.03 * yr;
    for i = 1:n
        stars = pStars(pval(i));
        if isempty(stars), continue; end
        if betas(i) >= 0
            text(i, up(i) + offset, stars, 'Color', 'r', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                 'FontWeight', 'bold');
        else
            text(i, lo(i) - offset, stars, 'Color', 'r', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                 'FontWeight', 'bold');
        end
    end

    set(gca, 'XTick', 1:n, 'XTickLabel', names, ...
             'XTickLabelRotation', 45, ...
             'TickLabelInterpreter', 'none');

    pad = 0.10 * yr;
    ylim([min(lo) - pad, max(up) + pad]);
    ylabel('\beta (fixed effect)');
    box off;
    hold off;
end

function s = pStars(p)
    if     p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = '';
    end
end
