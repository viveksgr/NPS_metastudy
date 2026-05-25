function [delta_matrix, subject_ids, subject_stats] = subject_condition_delta_matrix(X, subject, condition, posValue, baseValue)
% subject_condition_delta_matrix  Reduce trial-level matrix data to subject-wise deltas.
%
% [delta_matrix, subject_ids, subject_stats] = subject_condition_delta_matrix(X, subject, condition)
% computes, for each subject and each row of X:
%
%     mean(X(:, condition == 1), 2) - mean(X(:, condition == 0), 2)
%
% within that subject. Trials with condition -1 or any other value are ignored.
% Means ignore NaNs. This is intended for data shaped Nfeatures x Nt, such as
% Nvox x pooled trials fMRI matrices.
%
% Inputs:
%   X         : Nfeatures x Nt numeric matrix
%   subject   : 1 x Nt or Nt x 1 vector/cell array/string array of subject ids
%   condition : 1 x Nt or Nt x 1 numeric vector of condition codes
%   posValue  : optional condition value for the positive mean (default 1)
%   baseValue : optional condition value for the baseline mean (default 0)
%
% Outputs:
%   delta_matrix : Nfeatures x nSubjects matrix
%   subject_ids  : nSubjects x 1 subject identifiers, stable by first occurrence
%   subject_stats: table with per-subject condition trial counts

    if nargin < 4 || isempty(posValue)
        posValue = 1;
    end

    if nargin < 5 || isempty(baseValue)
        baseValue = 0;
    end

    validateattributes(X, {'numeric'}, {'2d'}, mfilename, 'X');
    validateattributes(condition, {'numeric'}, {'vector'}, mfilename, 'condition');

    subject = subject(:);
    condition = condition(:);

    nTrials = size(X, 2);
    if numel(subject) ~= nTrials || numel(condition) ~= nTrials
        error('X must have Nt columns, and subject/condition must both have Nt elements.');
    end

    subject_key = string(subject);
    [subject_key_ids, firstIdx] = unique(subject_key, 'stable');
    subject_ids = subject(firstIdx);
    [~, G] = ismember(subject_key, subject_key_ids);
    nSubjects = numel(subject_ids);
    nFeatures = size(X, 1);

    delta_matrix = nan(nFeatures, nSubjects);
    n_pos = zeros(nSubjects, 1);
    n_base = zeros(nSubjects, 1);

    for s = 1:nSubjects
        isSubject = G == s;
        posIdx = isSubject & condition == posValue;
        baseIdx = isSubject & condition == baseValue;

        n_pos(s) = sum(posIdx);
        n_base(s) = sum(baseIdx);

        pos_mean = mean(X(:, posIdx), 2, 'omitnan');
        base_mean = mean(X(:, baseIdx), 2, 'omitnan');
        delta_matrix(:, s) = pos_mean - base_mean;
    end

    subject_stats = table(subject_ids, n_pos, n_base, ...
        'VariableNames', {'subject', 'n_pos', 'n_base'});
end
