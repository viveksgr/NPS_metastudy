function [outMat, fieldNames] = collateContrastScalars(contrastCells, summaryFunc)
% COLLATECONTRASTSCALARS  Extract scalar fields from N contrast structs
%
%   [outMat, fieldNames] = collateContrastScalars(contrastCells, summaryFunc)
%
% INPUTS:
%   contrastCells – 1×N cell array of canlabcontrast (or similar) objects
%   summaryFunc   – function handle: s = summaryFunc(contrastCells{i});
%                   must return a scalar‐valued struct `s` with fields
%                   (one field named 'warnings' will be dropped if present)
%
% OUTPUTS:
%   outMat     – M×N array of the scalar values, where M = number of fields
%   fieldNames – M×1 cell array of the field names in `outMat` rows
%
% EXAMPLE:
%   % suppose you have
%   %   C = {c1, c2, …, cN};   % your canlabcontrast objects
%   %   f = @(c) myContrastSummary(c);
%   [A, names] = collateContrastScalars(C, f);
%   % then A(k,i) is the value of field names{k} for study i

N = numel(contrastCells);
if N==0
    outMat = []; fieldNames = {}; return;
end

% 1) Apply summaryFunc to each contrast, collect structs
structCells = cell(size(contrastCells));     % prealloc
for i = 1:N
    fname = contrastCells{i};
    fname = rescale(fname, 'l2norm_images');
    structCells{i} = summaryFunc( fname);
end


% 1) get and validate field names
allFields = fieldnames(structCells{1});
keep = ~strcmp(allFields, 'warnings');
fieldNames = allFields(keep);
M = numel(fieldNames);

% optional sanity: ensure every struct has same fields
for i=2:N
    if ~isequal(fieldnames(structCells{i}), allFields)
        error('Struct %d has different fields.', i);
    end
end

% 2) preallocate output
outMat = nan(M, N);

% 3) fill in
for i = 1:N
    S = structCells{i};
    for k = 1:M
        val = S.(fieldNames{k});
        if ~isnumeric(val) || ~isscalar(val)
            error('Field "%s" in struct %d is not a numeric scalar.', fieldNames{k}, i);
        end
        outMat(k,i) = val;
    end
end
