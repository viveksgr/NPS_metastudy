function [ DataCols] = harvest_canlab_funccols(rootdir, Idx)
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

for s = 1:S
    fprintf('Study %d/%d  (%s)\n', s, S, subdirs(s).name);

    % locate canlab_dataset*/canlab*.mat
    inner = dir(fullfile(rootdir, subdirs(s).name, 'fmri_data*'));
    if isempty(inner), error('No canlab_dataset* dir in %s', subdirs(s).name); end
    mfile = dir(fullfile(inner(1).folder, inner(1).name, 'contrast_data*.mat'));
    if isempty(mfile), error('No canlab*.mat in %s', inner(1).name); end

    load(fullfile(mfile(1).folder, mfile(1).name), 'DATA_OBJ_CON');   % loads DAT

    colIdx = Idx{s};
    ndats = length(colIdx);

    DataCols(end+1:end+ndats) = DATA_OBJ_CON(colIdx);        %#ok<*AGROW>
end



