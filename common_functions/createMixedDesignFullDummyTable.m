function [tbl, X, groupVec, varNames, info] = createMixedDesignFullDummyTable(studyCats, nsubjects, studyCovs, covNames)
% CREATEMIXEDDESIGNFULLDUMMYTABLE build subject-level table with one
% indicator column per study category.
%
% USAGE:
%   [tbl, X, groupVec, varNames, info] = createMixedDesignFullDummyTable(studyCats, nsubjects)
%   [tbl, X, groupVec, varNames, info] = createMixedDesignFullDummyTable(studyCats, nsubjects, studyCovs, covNames)
%
% INPUTS:
%   studyCats  - Sx1 categorical/cellstr/string array of study categories
%   nsubjects  - Sx1 numeric vector (# subjects per study)
%   studyCovs  - optional SxP matrix of study-level covariates
%   covNames   - optional 1xP cell array of names for studyCovs
%
% OUTPUTS:
%   tbl        - table with N = sum(nsubjects) rows and columns:
%                StudyID, StudyLabel, one dummy column per category,
%                and any expanded study-level covariates
%   X          - numeric design matrix matching tbl columns in varNames
%   groupVec   - integer vector mapping each row to study (1..S)
%   varNames   - names of the numeric columns in X, in order
%   info       - bookkeeping about levels and coding
%
% NOTES:
%   - For K category levels, this creates K dummy-coded columns.
%   - No reference level is omitted.
%   - Category columns are ordered alphabetically by class label.
%   - The generated fitlme formula uses -1 so each class beta estimates
%     that class mean directly rather than a difference from a reference.
%
% EXAMPLE:
%   [tbl, X, groupVec, varNames, info] = createMixedDesignFullDummyTable(Group_labels, nstud);
%   tbl.NPS = Y;
%   lme = fitlme(tbl, strrep(info.formula, 'Y', 'NPS'));

if nargin < 2
    error('Need studyCats and nsubjects.');
end

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

levels = sort(cellstr(unique(studyLabelVec)));

tbl = table();
tbl.StudyID = categorical(groupVec);
tbl.StudyLabel = categorical(cellstr(studyLabelVec), levels);

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

D = dummyvar(tbl.StudyLabel);
dummyColNames = cell(1, numel(levels));

for k = 1:numel(levels)
    rawName = sprintf('%s', levels{k});
    cleanName = matlab.lang.makeValidName(rawName);
    tbl.(cleanName) = D(:, k);
    dummyColNames{k} = cleanName;
end

X = D;
varNames = dummyColNames;

for q = 1:P
    X = [X, tbl.(covNames{q})]; %#ok<AGROW>
    varNames{end + 1} = covNames{q}; %#ok<AGROW>
end

info = struct();
info.coding = 'full-dummy/no-reference';
info.levels = levels;
info.referenceLevel = [];
info.referenceIndex = [];
info.dummyColNames = dummyColNames;
info.groupVec = groupVec;
info.varNames = varNames;
info.interceptMeaning = 'No fixed-effect intercept is included in the generated formula.';
info.betaMeaning = 'Each label beta is the estimated mean for that level.';

info.colIndex.dummies = 1:numel(dummyColNames);
info.colIndex.covariates = (numel(dummyColNames) + 1):(numel(dummyColNames) + P);

info.fullDummyCoding = D;
info.keepIdx = 1:numel(levels);

if isempty(dummyColNames)
    fixedPart = '-1';
else
    fixedPart = sprintf('-1 + %s', strjoin(dummyColNames, ' + '));
end

if ~isempty(covNames)
    fixedPart = sprintf('%s + %s', fixedPart, strjoin(covNames, ' + '));
end

info.formula = sprintf('Y ~ %s + (1|StudyID)', fixedPart);

end
