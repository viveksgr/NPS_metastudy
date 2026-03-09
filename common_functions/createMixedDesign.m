function [tbl, X, groupVec, varNames] = createMixedDesign(studyCats, nsubjects, studyCovs, covNames, varargin)
% CREATE MIXED-DESIGN build subject-level table + numeric design matrix
%
% USAGE:
%   [tbl, X, groupVec, varNames] = createMixedDesign(studyCats, nsubjects)
%   [tbl, X, groupVec, varNames] = createMixedDesign(studyCats, nsubjects, studyCovs, covNames)
%
% INPUTS:
%   studyCats  - Sx1 categorical/cell array (e.g. {'conditioning','suggestion',...})
%   nsubjects  - Sx1 numeric vector (# subjects per study)
%   studyCovs  - (optional) S x P numeric matrix of study-level continuous covariates
%   covNames   - (optional) 1xP cellstr names for studyCovs
%
% OUTPUTS:
%   tbl        - table with N = sum(nsubjects) rows: columns StudyID, (repeated) study-level predictors
%   X          - numeric design matrix (N x M) : [Intercept, categorical effectcols, studyCovs]
%   groupVec   - integer vector length N with study index for each row (1..S)
%   varNames   - cellstr names for columns of X (in order)
%
% NOTES:
%  - categorical predictor is effect-centered (columns sum to zero); it uses k-1 columns for k levels.
%  - If you pass no studyCovs, only categorical variable will be used (if provided).
%  - Use fitlme(tbl, 'Y ~ 1 + PredictorName + (1|StudyID)') to fit a mixed model.

if nargin < 2, error('Need studyCats and nsubjects.'); end
S = numel(nsubjects);
if numel(studyCats) ~= S, error('studyCats and nsubjects must match length S'); end

% default empty
if nargin < 3 || isempty(studyCovs), studyCovs = []; end
if nargin < 4 || isempty(covNames)
    if ~isempty(studyCovs)
        covNames = arrayfun(@(i) sprintf('cov%d',i), 1:size(studyCovs,2), 'Uni', false);
    else
        covNames = {};
    end
end

% 1) Expand rows to subjects
N = sum(nsubjects);
groupVec = zeros(N,1);
studyLabelVec = cell(N,1);
row = 1;
for s = 1:S
    ns = nsubjects(s);
    idx = row:(row+ns-1);
    groupVec(idx) = s;
    studyLabel = studyCats{s};
    if iscell(studyLabel), studyLabel = studyLabel{1}; end
    studyLabelVec(idx) = repmat({studyLabel}, ns, 1);
    row = row + ns;
end

% 2) Build table with StudyID and categorical predictor (study-level)
tbl = table();
tbl.StudyID = categorical(groupVec);  % grouping factor for fitlme
tbl.StudyLabel = categorical(studyLabelVec); % study-level category (e.g., 'conditioning' etc.)

% 3) Expand studyCovs to subject-level
P = size(studyCovs,2);
for p = 1:P
    colvals = nan(N,1);
    row = 1;
    for s = 1:S
        ns = nsubjects(s);
        colvals(row:row+ns-1) = studyCovs(s,p);
        row = row + ns;
    end
    tbl.(covNames{p}) = colvals;
end

% 4) Make numeric design matrix X
% Intercept
X = ones(N,1);
varNames = {'Intercept'};

% 4a) categorical study label -> dummy/contrast columns (effect-centered)
levels = categories(tbl.StudyLabel);
K = numel(levels);
if K > 1
    D = dummyvar(tbl.StudyLabel);      % N x K  (each column = 1 for that level)
    % drop last column to avoid linear dependence and center columns to sum to 0
    Dc = D(:,1:K-1);
    Dc = Dc - mean(Dc,1);  % effect-center so each column sums to zero
    X = [X Dc];
    for k = 1:(K-1), varNames{end+1} = sprintf('Label_%s', levels{k}); end
end

% 4b) continuous study covariates (already expanded in tbl)
for p = 1:P
    X = [X tbl.(covNames{p})];
    varNames{end+1} = covNames{p};
end

end
