function [p_thresh, p_correct] = fdr_benjhoc(p_vals, q)
% FDR_BENJHOC   Benjamini–Hochberg FDR correction
%
%   [p_thresh, p_correct] = fdr_benjhoc(p_vals, q)
%
% INPUTS:
%   p_vals   – [N×1] vector of uncorrected p-values
%   q        – desired FDR rate (scalar, default 0.05)
%
% OUTPUTS:
%   p_thresh – the largest uncorrected p such that
%              p_vals <= p_thresh ⇒ reject H0 under FDR<q
%   p_correct– [N×1] vector of FDR-adjusted p-values (≤1)
%
% Example:
%   [pth, p_adj] = fdr_benjhoc(pvals, 0.1);

    if nargin<2 || isempty(q), q = 0.05; end
    p = p_vals(:);
    N = numel(p);
    if N==0
        p_thresh = []; p_correct = [];
        return
    end

    % Sort p-values ascending
    [p_sorted, sort_idx] = sort(p);
    rank = (1:N)';

    % Find threshold index
    crit = (rank / N) * q;
    below = find(p_sorted <= crit, 1, 'last');
    if isempty(below)
        p_thresh = 0;
    else
        p_thresh = p_sorted(below);
    end

    % Compute raw adjusted p: p_sorted .* N ./ rank
    raw_adj = p_sorted * N ./ rank;

    % Enforce monotonicity:  
    % adjusted p at i is min(raw_adj(i:end))
    % and cap at 1
    p_adj_sorted = min( cummin(raw_adj(end:-1:1)), 1 );
    p_adj_sorted = p_adj_sorted(end:-1:1);

    % Unscramble to original order
    p_correct = nan(N,1);
    p_correct(sort_idx) = p_adj_sorted;
end
