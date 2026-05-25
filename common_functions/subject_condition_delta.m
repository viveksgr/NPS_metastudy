function [delta, subject_ids, subject_stats] = subject_condition_delta(y, subject, condition, posValue, baseValue)
% subject_condition_delta  Reduce trial-level values to subject-wise condition deltas.
%
% [delta, subject_ids, subject_stats] = subject_condition_delta(y, subject, condition)
% computes, for each subject:
%
%     mean(y(condition == 1)) - mean(y(condition == 0))
%
% Trials with condition -1 or any other value are ignored. Means ignore NaNs.
%
% Inputs:
%   y         : 1 x Nt or Nt x 1 numeric vector of trial-level values
%   subject   : 1 x Nt or Nt x 1 vector/cell array/string array of subject ids
%   condition : 1 x Nt or Nt x 1 numeric vector of condition codes
%   posValue  : optional condition value for the positive mean (default 1)
%   baseValue : optional condition value for the baseline mean (default 0)
%
% Outputs:
%   delta        : 1 x nSubjects vector, in subject_ids order
%   subject_ids  : nSubjects x 1 subject identifiers, stable by first occurrence
%   subject_stats: table with per-subject counts, means, and delta

    if nargin < 4 || isempty(posValue)
        posValue = 1;
    end

    if nargin < 5 || isempty(baseValue)
        baseValue = 0;
    end

    validateattributes(y, {'numeric'}, {'vector'}, mfilename, 'y');
    validateattributes(condition, {'numeric'}, {'vector'}, mfilename, 'condition');

    y = y(:);
    condition = condition(:);
    subject = subject(:);

    nTrials = numel(y);
    if numel(subject) ~= nTrials || numel(condition) ~= nTrials
        error('y, subject, and condition must have the same number of trials.');
    end

    subject_key = string(subject);
    [subject_key_ids, firstIdx] = unique(subject_key, 'stable');
    subject_ids = subject(firstIdx);
    [~, G] = ismember(subject_key, subject_key_ids);
    nSubjects = numel(subject_ids);

    pos_mean = nan(nSubjects, 1);
    base_mean = nan(nSubjects, 1);
    n_pos = zeros(nSubjects, 1);
    n_base = zeros(nSubjects, 1);

    for s = 1:nSubjects
        isSubject = G == s;
        posIdx = isSubject & condition == posValue;
        baseIdx = isSubject & condition == baseValue;

        n_pos(s) = sum(posIdx);
        n_base(s) = sum(baseIdx);
        pos_mean(s) = mean(y(posIdx), 'omitnan');
        base_mean(s) = mean(y(baseIdx), 'omitnan');
    end

    delta_col = pos_mean - base_mean;
    delta = delta_col.';

    subject_stats = table(subject_ids, n_pos, n_base, pos_mean, base_mean, delta_col, ...
        'VariableNames', {'subject', 'n_pos', 'n_base', 'pos_mean', 'base_mean', 'delta'});
end
