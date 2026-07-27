function [tbl, X, groupVec, varNames, info] = createMixedDesignDummyTable(studyCats, nsubjects, studyCovs, covNames, varargin)
% CREATEMIXEDDESIGNDUMMYTABLE build subject-level table with explicit
% dummy/reference-coded regressors.
%
% USAGE:
%   [tbl, X, groupVec, varNames, info] = createMixedDesignDummyTable(studyCats, nsubjects)
%   [tbl, X, groupVec, varNames, info] = createMixedDesignDummyTable(studyCats, nsubjects, studyCovs, covNames)
%   [tbl, X, groupVec, varNames, info] = createMixedDesignDummyTable(..., 'referenceLevel', 'Cognitive')
%
% INPUTS:
%   studyCats  - Sx1 categorical/cellstr/string array of study categories
%   nsubjects  - Sx1 numeric vector (# subjects per study)
%   studyCovs  - optional SxP matrix of study-level covariates
%   covNames   - optional 1xP cell array of names for studyCovs
%
% OPTIONAL NAME/VALUE:
%   'referenceLevel' - category label to omit from dummy-coded columns.
%                      Default is the last category level.
%
% OUTPUTS:
%   tbl        - table with N = sum(nsubjects) rows and columns:
%                StudyID, StudyLabel, dummy-coded regressors,
%                and any expanded study-level covariates
%   X          - numeric design matrix matching tbl columns in varNames
%   groupVec   - integer vector mapping each row to study (1..S)
%   varNames   - names of the numeric columns in X, in order
%   info       - bookkeeping about levels and coding
%
% NOTES:
%   - For K category levels, this creates K-1 dummy-coded columns.
%   - The reference level is omitted.
%   - With an intercept in the model:
%       Intercept = estimated mean of reference level
%       Beta_k    = estimated difference: level k - reference level
%
% EXAMPLE:
%   [tbl, X, groupVec, varNames, info] = createMixedDesignDummyTable(Group_labels, nstud, ...
%       [], [], 'referenceLevel', 'Cognitive');
%   tbl.NPS = Y;
%   lme = fitlme(tbl, strrep(info.formula, 'Y', 'NPS'));

if nargin < 2
    error('Need studyCats and nsubjects.');
end

p = inputParser;
addParameter(p, 'referenceLevel', [], @(x) ischar(x) || isstring(x) || iscellstr(x) || isempty(x));
parse(p, varargin{:});
referenceLevel = p.Results.referenceLevel;

S = numel(nsubjects);
if numel(studyCats) ~= S
    error('studyCats and nsubjects must match length S.');
end

if nargin < 3 || isempty(studyCovs)
    studyCovs = [];
end

if nargin < 4 || isempty(covNames)
    if ~isempty(studyCovs)
        covNames = arrayfun(@(i) sprintf('cov%d', i), 1:size(studyCovs, 2), 'UniformOutput', false);
    else
        covNames = {};
    end
end

if ~isempty(studyCovs) && size(studyCovs, 1) ~= S
    error('studyCovs must have one row per study.');
end

if ~isempty(studyCovs) && numel(covNames) ~= size(studyCovs, 2)
    error('covNames must match the number of columns in studyCovs.');
end

nsubjects = nsubjects(:);
N = sum(nsubjects);

groupVec = zeros(N, 1);
studyLabelVec = strings(N, 1);

row = 1;
for s = 1:S
    ns = nsubjects(s);
    idx = row:(row + ns - 1);
    groupVec(idx) = s;

    label = studyCats(s);
    if iscell(label)
        label = label{1};
    end
    label = string(label);
    studyLabelVec(idx) = repmat(label, ns, 1);

    row = row + ns;
end

tbl = table();
tbl.StudyID = categorical(groupVec);
tbl.StudyLabel = categorical(cellstr(studyLabelVec));

% Expand study-level covariates to subject-level.
P = size(studyCovs, 2);
for q = 1:P
    colvals = nan(N, 1);
    row = 1;
    for s = 1:S
        ns = nsubjects(s);
        colvals(row:(row + ns - 1)) = studyCovs(s, q);
        row = row + ns;
    end
    tbl.(covNames{q}) = colvals;
end

X = ones(N, 1);
varNames = {'Intercept'};

levels = categories(tbl.StudyLabel);
K = numel(levels);
dummyColNames = {};

if K > 1
    if isempty(referenceLevel)
        refIdx = K;
        referenceLevel = levels{refIdx};
    else
        referenceLevel = char(string(referenceLevel));
        refIdx = find(strcmp(levels, referenceLevel), 1);
        if isempty(refIdx)
            error('referenceLevel "%s" was not found in studyCats.', referenceLevel);
        end
    end

    D = dummyvar(tbl.StudyLabel);      % N x K, pure dummy coding
    keepIdx = setdiff(1:K, refIdx, 'stable');
    DReduced = D(:, keepIdx);          % omit reference category

    for k = 1:numel(keepIdx)
        rawName = sprintf('%s', levels{keepIdx(k)});
        cleanName = matlab.lang.makeValidName(rawName);
        tbl.(cleanName) = DReduced(:, k);
        dummyColNames{end + 1} = cleanName; %#ok<AGROW>
    end

    X = [X, DReduced];
    varNames = [varNames, dummyColNames];

else
    refIdx = 1;
    referenceLevel = levels{1};
    D = ones(N, 1);
    DReduced = zeros(N, 0);
end

for q = 1:P
    X = [X, tbl.(covNames{q})];
    varNames{end + 1} = covNames{q}; %#ok<AGROW>
end

info = struct();
info.coding = 'dummy/reference';
info.levels = levels;
info.referenceLevel = referenceLevel;
info.referenceIndex = refIdx;
info.dummyColNames = dummyColNames;
info.groupVec = groupVec;
info.varNames = varNames;
info.interceptMeaning = 'Intercept is the estimated mean of the reference level.';
info.betaMeaning = 'Each label beta is the estimated difference between that level and the reference level.';

info.colIndex.intercept = 1;
info.colIndex.dummies = 2:(1 + numel(dummyColNames));
info.colIndex.covariates = (2 + numel(dummyColNames)):(1 + numel(dummyColNames) + P);

info.fullDummyCoding = D;
info.reducedDummyCoding = DReduced;
info.keepIdx = keepIdx;

if isempty(dummyColNames)
    fixedPart = '1';
else
    fixedPart = sprintf('1 + %s', strjoin(dummyColNames, ' + '));
end

if ~isempty(covNames)
    fixedPart = sprintf('%s + %s', fixedPart, strjoin(covNames, ' + '));
end

info.formula = sprintf('Y ~ %s + (1|StudyID)', fixedPart);

end