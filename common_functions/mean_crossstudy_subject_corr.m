function [mean_r, se_r, mean_z, se_z, nPairs] = mean_crossstudy_subject_corr(studies, varargin)
% MEAN_CROSSSTUDY_SUBJECT_CORR  Mean subject–subject correlation of each study
%
% [mean_r, se_r, mean_z, se_z, nPairs] = mean_crossstudy_subject_corr(studies, ...)
%
% Input:
%   studies  - 1xS cell array. Each cell is nVox x nSubs numeric matrix
%              OR an object/struct with .dat (nVox x nSubs).
%
% Name/value pairs:
%   'useFisher' (true)  - average in Fisher z-space (recommended)
%   'verbose'   (false) - show progress
%   'blockSize' (0)     - if >0, compute correlations in blocks of this many columns to save memory
%
% Outputs:
%   mean_r  - Sx1 mean correlation (back-transformed if useFisher true)
%   se_r    - Sx1 approximate SE of mean_r (delta method from z-space)
%   mean_z  - Sx1 mean Fisher z used (atanh of r)
%   se_z    - Sx1 SE in z-space
%   nPairs  - Sx1 number of subject pairs used per study (sum across other studies)
%
% Example:
%   [mr, se] = mean_crossstudy_subject_corr(myStudiesCell, 'blockSize', 500);

% parse inputs
p = inputParser;
addParameter(p,'useFisher', true, @islogical);
addParameter(p,'verbose', false, @islogical);
addParameter(p,'blockSize', 0, @(x) isnumeric(x) && isscalar(x) && x>=0);
parse(p, varargin{:});
useFisher = p.Results.useFisher;
doVerbose = p.Results.verbose;
blockSize = round(p.Results.blockSize);

S = numel(studies);
if S < 2
    error('Need at least two studies in the cell array.');
end

% extract data matrices & check consistent voxel count
Xcell = cell(1,S);
nV = NaN;
for s = 1:S
    cur = studies{s};
    if isstruct(cur) || isobject(cur)
        if isfield(cur,'dat') || isprop(cur,'dat')
            X = double(cur.dat);
        else
            error('Study %d is object/struct but has no .dat field.', s);
        end
    else
        X = double(cur);
    end
    if isempty(X)
        error('Study %d has empty data.', s);
    end
    [nv, ns] = size(X);
    if isnan(nV)
        nV = nv;
    elseif nv ~= nV
        error('Voxel count mismatch: study %d has %d voxels but expected %d.', s, nv, nV);
    end
    Xcell{s} = X;
end

% prepare outputs
mean_r = nan(S,1);
se_r   = nan(S,1);
mean_z = nan(S,1);
se_z   = nan(S,1);
nPairs = zeros(S,1);

% main loop: for each study i, correlate its subjects with all other subjects
for i = 1:S
    if doVerbose, fprintf('Study %d / %d\n', i, S); end
    Xi = Xcell{i};            % nV x ni
    [~, ni] = size(Xi);
    if ni < 1
        if doVerbose, fprintf('  study %d has zero subjects, skipping\n', i); end
        continue;
    end

    % build list (or stream) of other-study matrices
    otherIdx = setdiff(1:S, i);
    % total columns in other studies
    totalOthers = sum(cellfun(@(c) size(c,2), Xcell(otherIdx)));

    % if blockSize==0 (no blocking) we'll concatenate all others then corr in one call
    % else process in blocks across other columns
    zvals_accum = []; % will store Fisher z values (if useFisher) or r values (if not)
    rcount = 0;

    if blockSize <= 0
        % Concatenate all others
        Xothers = [];
        for j = otherIdx
            Xothers = [Xothers, Xcell{j}]; %#ok<AGROW>
        end
        if isempty(Xothers)
            if doVerbose, fprintf('  no other subjects for study %d\n', i); end
            mean_r(i) = NaN; se_r(i) = NaN; mean_z(i) = NaN; se_z(i) = NaN; nPairs(i)=0;
            continue;
        end
        % compute correlations between columns of Xi and columns of Xothers
        % corr(Xi, Xothers, 'Rows','pairwise') gives ni x Nothers
        R = corr(Xi, Xothers, 'Rows', 'pairwise');
        rvals = R(:);
        valid = ~isnan(rvals);
        rvals = rvals(valid);
        rcount = numel(rvals);
        if rcount == 0
            if doVerbose, fprintf('  no valid subject–subject pairs for study %d (all NaN)\n', i); end
            mean_r(i) = NaN; se_r(i) = NaN; mean_z(i) = NaN; se_z(i) = NaN; nPairs(i)=0;
            continue;
        end
        % clamp extremes
        clamp = 1 - 1e-6;
        rvals(rvals >= clamp) = clamp;
        rvals(rvals <= -clamp) = -clamp;

        if useFisher
            zvals = atanh(rvals);
            mean_z(i) = mean(zvals, 'omitnan');
            se_z(i)   = std(zvals, 0, 'omitnan') / sqrt(numel(zvals));
            mean_r(i) = tanh(mean_z(i));
            se_r(i)   = (1 - mean_r(i).^2) * se_z(i); % delta-method
        else
            mean_r(i) = mean(rvals, 'omitnan');
            se_r(i)   = std(rvals, 0, 'omitnan') / sqrt(numel(rvals));
            mean_z(i) = atanh(mean_r(i));
            se_z(i)   = NaN;
        end
        nPairs(i) = rcount;

    else
        % BLOCKING mode: iterate through other-study matrices in blocks of columns
        % Build a list of all other columns as sequential concatenation of other studies
        % but avoid building full concatenation at once.
        % Create an index vector of column counts per other study
        otherColsMat = cell2mat(cellfun(@(c) size(c,2), Xcell(otherIdx), 'UniformOutput', false));
        % We'll iterate through each other study sequentially, but within each other study we may block columns
        for j = otherIdx
            Xj = Xcell{j};  % nV x nj
            nj = size(Xj,2);
            if nj == 0, continue; end
            if blockSize >= nj
                % compute full corr Xi x Xj
                Rblk = corr(Xi, Xj, 'Rows','pairwise');
                rvals = Rblk(:);
                valid = ~isnan(rvals);
                rvals = rvals(valid);
            else
                % block through columns of Xj
                rvals = [];
                startc = 1;
                while startc <= nj
                    endc = min(nj, startc + blockSize - 1);
                    Rblk = corr(Xi, Xj(:, startc:endc), 'Rows','pairwise');
                    rr = Rblk(:);
                    rr = rr(~isnan(rr));
                    rvals = [rvals; rr]; %#ok<AGROW>
                    startc = endc + 1;
                end
            end
            if isempty(rvals), continue; end
            % clamp and transform incremental accumulation
            clamp = 1 - 1e-6;
            rvals(rvals >= clamp) = clamp;
            rvals(rvals <= -clamp) = -clamp;
            if useFisher
                zvals_accum = [zvals_accum; atanh(rvals)]; %#ok<AGROW>
            else
                zvals_accum = [zvals_accum; rvals]; %#ok<AGROW>
            end
            rcount = rcount + numel(rvals);
        end

        if rcount == 0
            if doVerbose, fprintf('  no valid pairs for study %d\n', i); end
            mean_r(i)=NaN; se_r(i)=NaN; mean_z(i)=NaN; se_z(i)=NaN; nPairs(i)=0;
            continue;
        end

        if useFisher
            mean_z(i) = mean(zvals_accum, 'omitnan');
            se_z(i)   = std(zvals_accum, 0, 'omitnan') ./ sqrt(numel(zvals_accum));
            mean_r(i) = tanh(mean_z(i));
            se_r(i)   = (1 - mean_r(i).^2) * se_z(i);
        else
            mean_r(i) = mean(zvals_accum, 'omitnan');
            se_r(i)   = std(zvals_accum, 0, 'omitnan') ./ sqrt(numel(zvals_accum));
            mean_z(i) = atanh(mean_r(i));
            se_z(i)   = NaN;
        end
        nPairs(i) = rcount;
    end

    if doVerbose
        fprintf('  study %d: pairs=%d, mean_r=%.4f, se_r=%.4f\n', i, nPairs(i), mean_r(i), se_r(i));
    end
end

end
