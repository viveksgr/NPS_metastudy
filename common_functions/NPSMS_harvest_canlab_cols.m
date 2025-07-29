function [BigDat, ColNames] = NPSMS_harvest_canlab_cols(rootdir, Idx, roughMaxRows)
% HARVEST_CANLAB_COLS  Concatenate selected Subj_Level columns from many studies
%
%  [BigDat, ColNames] = harvest_canlab_cols(rootdir, Idx, roughMaxRows)
%
%  INPUTS
%    rootdir       – directory whose immediate sub-folders are the studies
%    Idx           – 1×S cell; Idx{s} is a vector of column indices to grab
%                    from DAT.Subj_Level for study s
%    roughMaxRows  – (optional) rough upper bound of subjects; default 100
%
%  OUTPUTS
%    BigDat   – [Nmax × Nt] matrix, NaN-padded, where Nt = Σ numel(Idx{s})
%    ColNames – 1×Nt cell array of contrast names in the same column order
%
%  EXAMPLE
%    [BD, names] = harvest_canlab_cols('C:\MyStudies', Idx, 120);

if nargin < 3 || isempty(roughMaxRows), roughMaxRows = 100; end

% -------------------------------------------------------------------------
% 0.  Enumerate studies (sub-folders of rootdir) in alphabetical order
% -------------------------------------------------------------------------
d = dir(rootdir);
subdirs = d([d.isdir] & ~startsWith({d.name}, '.'));
subdirs = sortrows(struct2table(subdirs), 'name');
subdirs = table2struct(subdirs);
S = numel(subdirs);

if S ~= numel(Idx)
    error('Idx has %d entries, but %d subdirectories found.', numel(Idx), S);
end

% -------------------------------------------------------------------------
% 1.  Loop through studies, extract requested columns
% -------------------------------------------------------------------------
DataCols  = {};
NameCols  = {};
maxRows   = 0;

for s = 1:S
    if s==3
        'beep'
    end
    fprintf('Study %d/%d  (%s)\n', s, S, subdirs(s).name);

    % locate canlab_dataset*/canlab*.mat
    inner = dir(fullfile(rootdir, subdirs(s).name, 'canlab_dataset*'));
    if isempty(inner), error('No canlab_dataset* dir in %s', subdirs(s).name); end
    mfile = dir(fullfile(inner(1).folder, inner(1).name, 'canlab*.mat'));
    if isempty(mfile), error('No canlab*.mat in %s', inner(1).name); end

    load(fullfile(mfile(1).folder, mfile(1).name), 'DAT');   % loads DAT

    colIdx = Idx{s};
    if max(colIdx) > size(DAT.Subj_Level.data,2)
        error('Idx{%d} requests column %d but study has only %d columns.', ...
                s, max(colIdx), size(DAT.Subj_Level.data,2));
    end

    if colIdx>0
    DataCols{end+1} = DAT.Subj_Level.data(:, colIdx);        %#ok<*AGROW>
    NameCols{end+1} = DAT.Subj_Level.names(colIdx);
    maxRows = max(maxRows, size(DataCols{end},1));
    end
end

% allow for larger studies than roughMaxRows
maxRows = max(maxRows, roughMaxRows);

% -------------------------------------------------------------------------
% 2.  Assemble BigDat (NaN-padded) and ColNames
% -------------------------------------------------------------------------
totalCols = sum(cellfun(@(x) size(x,2), DataCols));
BigDat    = nan(maxRows, totalCols);
ColNames  = cell(1, totalCols);

colStart = 1;
for s = 1:S
    this = DataCols{s};
    r    = size(this,1);
    c    = size(this,2);
    BigDat(1:r, colStart:colStart+c-1) = this;
    ColNames(colStart:colStart+c-1)    = NameCols{s};
    colStart = colStart + c;
end

% -------------------------------------------------------------------------
% 3.  Trim trailing all-NaN rows (optional)
% -------------------------------------------------------------------------
lastRow = find(any(~isnan(BigDat),2), 1, 'last');
if ~isempty(lastRow)
    BigDat = BigDat(1:lastRow, :);
else
    warning('All rows are NaN — check your inputs.');
end
end
