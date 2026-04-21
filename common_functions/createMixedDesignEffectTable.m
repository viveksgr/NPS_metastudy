function [tbl, X, groupVec, varNames, info] = createMixedDesignEffectTable(studyCats, nsubjects, studyCovs, covNames, varargin)
% CREATEMIXEDDESIGNEFFECTTABLE build subject-level table with explicit
% effect-coded regressors so the intercept is the grand mean.
%
% USAGE:
%   [tbl, X, groupVec, varNames, info] = createMixedDesignEffectTable(studyCats, nsubjects)
%   [tbl, X, groupVec, varNames, info] = createMixedDesignEffectTable(studyCats, nsubjects, studyCovs, covNames)
%   [tbl, X, groupVec, varNames, info] = createMixedDesignEffectTable(..., 'referenceLevel', 'cognitive')
%
% INPUTS:
%   studyCats  - Sx1 categorical/cellstr/string array of study categories
%   nsubjects  - Sx1 numeric vector (# subjects per study)
%   studyCovs  - optional SxP matrix of study-level covariates
%   covNames   - optional 1xP cell array of names for studyCovs
%
% OPTIONAL NAME/VALUE:
%   'referenceLevel' - category label to omit from the explicit effect-coded
%                      columns. Default is the last category level.
%
% OUTPUTS:
%   tbl        - table with N = sum(nsubjects) rows and columns:
%                StudyID, StudyLabel, explicit effect-coded regressors,
%                and any expanded study-level covariates
%   X          - numeric design matrix matching tbl columns in varNames
%   groupVec   - integer vector mapping each row to study (1..S)
%   varNames   - names of the numeric columns in X, in order
%   info       - bookkeeping about levels and coding
%
% NOTES:
%   - For K category levels, this creates K-1 effect-coded columns.
%   - The omitted Kth level is implied by the negative sum of the others.
%   - With an intercept in the model, this parameterization makes the
%     intercept correspond to the row-weighted grand mean implied by the
%     expanded table.
%
% EXAMPLE:
%   [tbl, X, groupVec, varNames, info] = createMixedDesignEffectTable(Group_labels, nstud);
%   tbl.Y = Y;
%   lme = fitlme(tbl, info.formula);

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

% Expand study-level covariates to subject-level first so the output table
% contains all regressors needed by fitlme.
P = size(studyCovs, 2);
for p = 1:P
    colvals = nan(N, 1);
    row = 1;
    for s = 1:S
        ns = nsubjects(s);
        colvals(row:(row + ns - 1)) = studyCovs(s, p);
        row = row + ns;
    end
    tbl.(covNames{p}) = colvals;
end

X = ones(N, 1);
varNames = {'Intercept'};

levels = categories(tbl.StudyLabel);
K = numel(levels);
effectColNames = {};

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

    D = dummyvar(tbl.StudyLabel);      % N x K
    pLevel = mean(D, 1);               % observed level proportions
    Dc = D - pLevel;                   % full centered coding
    keepIdx = setdiff(1:K, refIdx, 'stable');
    DcReduced = Dc(:, keepIdx);        % drop chosen reference column for full rank

    for k = 1:numel(keepIdx)
        rawName = sprintf('Label_%s', levels{keepIdx(k)});
        cleanName = matlab.lang.makeValidName(rawName);
        tbl.(cleanName) = DcReduced(:, k);
        effectColNames{end + 1} = cleanName; %#ok<AGROW>
    end

    X = [X, DcReduced];
    varNames = [varNames, effectColNames];
else
    pLevel = 1;
    Dc = ones(N, 1) - 1;
    DcReduced = zeros(N, 0);
    refIdx = 1;
    referenceLevel = levels{1};
end

for p = 1:P
    X = [X, tbl.(covNames{p})];
    varNames{end + 1} = covNames{p}; %#ok<AGROW>
end

info = struct();
info.levels = levels;
info.referenceLevel = referenceLevel;
info.referenceIndex = refIdx;
info.effectColNames = effectColNames;
info.levelProportions = pLevel;
info.groupVec = groupVec;
info.varNames = varNames;
info.interceptMeaning = ['Intercept is the row-weighted grand mean under ' ...
    'this expanded-table coding.'];
info.colIndex.intercept = 1;
info.colIndex.effects = 2:(1 + numel(effectColNames));
info.colIndex.covariates = (2 + numel(effectColNames)):(1 + numel(effectColNames) + P);
info.fullCenteredCoding = Dc;
info.reducedCenteredCoding = DcReduced;

if isempty(effectColNames)
    info.formula = 'Y ~ 1 + (1|StudyID)';
else
    info.formula = sprintf('Y ~ 1 + %s + (1|StudyID)', strjoin(effectColNames, ' + '));
end

end
