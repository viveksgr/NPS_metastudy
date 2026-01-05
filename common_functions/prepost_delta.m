function [delta, subj_ids, subj_means] = prepost_delta(T, nameCol, stageCol, avgCol, deltaMode)
% prepost_delta  Compute subject-wise pre/post averages and delta.
%
% Inputs:
%   T         : table with columns for subject id, stage ('pre'/'post'), and values
%   nameCol   : string/char, column name for subject id (e.g., "name")
%   stageCol  : string/char, column name for stage (e.g., "stage")
%   avgCol    : string/char, column name for value (e.g., "avg")
%   deltaMode : (optional) "post-pre" (default) or "pre-post"
%
% Outputs:
%   delta     : num_subjects x 1 vector of (post - pre) (or reversed)
%   subj_ids  : subject identifiers in the same order as delta
%   subj_means: table with columns: subject, pre_mean, post_mean, delta

    if nargin < 5 || isempty(deltaMode)
        deltaMode = "post-pre";
    end

    % Extract columns
    subj  = T.(nameCol);
    stage = string(T.(stageCol));
    y     = T.(avgCol);

    % Normalize stage labels
    stage = lower(strtrim(stage));

    % Validate
    if ~all(ismember(stage, ["pre","post"]))
        bad = unique(stage(~ismember(stage, ["pre","post"])));
        error("Unexpected stage labels found: %s", strjoin(bad, ", "));
    end

    % Group by subject
    [G, subj_ids] = findgroups(subj);

    % Compute subject-wise means for each stage
    pre_mean  = splitapply(@(vals, st) mean(vals(st=="pre"),  'omitnan'),  y, stage, G);
    post_mean = splitapply(@(vals, st) mean(vals(st=="post"), 'omitnan'),  y, stage, G);

    % Delta
    switch string(deltaMode)
        case "post-pre"
            delta = post_mean - pre_mean;
        case "pre-post"
            delta = pre_mean - post_mean;
        otherwise
            error('deltaMode must be "post-pre" or "pre-post".');
    end

    % Package means
    subj_means = table(subj_ids, pre_mean, post_mean, delta, ...
        'VariableNames', {'subject','pre_mean','post_mean','delta'});
end
