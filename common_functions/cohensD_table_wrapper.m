function outTbl = cohensD_table_wrapper(data, varargin)
% COHENSD_TABLE_WRAPPER Compute Cohen's d (one-sample) for columns with varying lengths
%
% outTbl = cohensD_table_wrapper(data, 'nBoot',2000,'alpha',0.05,'verbose',false)
%
% Inputs:
%  - data : numeric matrix (rows padded with NaN), or cell array of numeric vectors,
%           or table (numeric columns with NaN padding, or columns that are cells holding vectors).
% Options:
%  - 'nBoot'  : number of bootstrap samples (default 2000)
%  - 'alpha'  : significance alpha for CI (default 0.05)
%  - 'verbose': true/false (default false)
%
% Output:
%  - outTbl : table with rows = columns/variables and fields:
%       Name    n   d    g    SE    CI_anal_lo CI_anal_hi   CI_boot_lo CI_boot_hi
%
% Example:
%   C = {randn(12,1)*0.5 + 0.3, randn(8,1)*0.8 + 0.1, randn(20,1)*0.2};
%   T = cohensD_table_wrapper(C,'nBoot',3000,'verbose',true)

% parse options
p = inputParser;
addParameter(p,'nBoot',2000,@(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p,'alpha',0.05,@(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
addParameter(p,'verbose',false,@islogical);
parse(p,varargin{:});
nBoot = p.Results.nBoot;
alpha = p.Results.alpha;
verbose = p.Results.verbose;

% Convert input into cell array of numeric vectors and names
[varCells, varNames] = normalize_input(data);

nVar = numel(varCells);

% Preallocate results
Ns = zeros(nVar,1);
ds = nan(nVar,1);
gs = nan(nVar,1);
SEs = nan(nVar,1);
CIa_lo = nan(nVar,1);
CIa_hi = nan(nVar,1);
CIB_lo = nan(nVar,1);
CIB_hi = nan(nVar,1);

% main loop: one-sample d per column (vs 0)
for k = 1:nVar
    x = varCells{k}(:);
    x = x(~isnan(x));  % drop NaNs
    n = numel(x);
    Ns(k) = n;
    if n < 2
        if verbose
            warning('Variable %s has n<2 (n=%d); results set to NaN', varNames{k}, n);
        end
        continue
    end

    % sample mean and sd
    m = mean(x);
    s = std(x,0);   % sample sd (n-1)
    d = m / s;

    % Hedges' small sample correction
    df = n - 1;
    J = 1 - (3 / (4*df - 1));
    g = d * J;

    % analytic SE (approx)
    SE = sqrt( (1/n) + (d.^2 ./ (2*n)) );

    % analytic CI
    z = norminv(1 - alpha/2);
    CI_anal = [d - z*SE, d + z*SE];

    % bootstrap percentile CI for d
    bootd = nan(nBoot,1);
    for b = 1:nBoot
        idx = randsample(n, n, true);
        xb = x(idx);
        mb = mean(xb);
        sb = std(xb,0);
        if sb == 0
            bootd(b) = NaN;
        else
            bootd(b) = mb / sb;
        end
    end
    bootd = bootd(~isnan(bootd));
    if isempty(bootd)
        CIB = [NaN NaN];
    else
        lo = 100*(alpha/2); hi = 100*(1-alpha/2);
        CIB = prctile(bootd, [lo hi]);
    end

    % store (report Cohen's d; Hedges' g returned too)
    ds(k) = d;
    gs(k) = g;
    SEs(k) = SE;
    CIa_lo(k) = CI_anal(1);
    CIa_hi(k) = CI_anal(2);
    CIB_lo(k) = CIB(1);
    CIB_hi(k) = CIB(2);

    if verbose
        fprintf('Var %s: n=%d, d=%.3f, g=%.3f, SE=%.3f, CI_anal=[%.3f %.3f], CI_boot=[%.3f %.3f]\n', ...
            varNames{k}, n, d, g, SE, CI_anal(1), CI_anal(2), CIB(1), CIB(2));
    end
end

% Build output table
outTbl = table(varNames(:), Ns, ds, gs, SEs, CIa_lo, CIa_hi, CIB_lo, CIB_hi, ...
    'VariableNames', {'Name','n','d','g','SE','CI_anal_lo','CI_anal_hi','CI_boot_lo','CI_boot_hi'});

end

%%%%%%%%%%% helper: normalize_input %%%%%%%%%%%
function [cellsOut, namesOut] = normalize_input(data)
% returns cell array of numeric vectors and names

if istable(data)
    T = data;
    names = T.Properties.VariableNames;
    % check if each column is a cell vector (e.g. single-row table with each cell containing vector)
    isCellCol = false(1,numel(names));
    for i = 1:numel(names)
        col = T.(names{i});
        if iscell(col) && numel(col)==1 && isnumeric(col{1})
            isCellCol(i) = true;
        end
    end
    if all(isCellCol)
        % table like: 1xN where each var contains a numeric vector in its only cell
        cellsOut = cell(1,numel(names));
        for i = 1:numel(names)
            cellsOut{i} = T.(names{i}){1};
        end
        namesOut = names;
        return
    else
        % assume numeric columns (possibly padded with NaN)
        cellsOut = cell(1,numel(names));
        for i = 1:numel(names)
            col = T.(names{i});
            if ~isnumeric(col)
                error('Table column %s must be numeric vector or a cell containing a numeric vector.', names{i});
            end
            cellsOut{i} = col(:);
        end
        namesOut = names;
        return
    end

elseif iscell(data)
    % cell array of numeric vectors
    nc = numel(data);
    cellsOut = cell(1,nc);
    namesOut = cell(1,nc);
    for i = 1:nc
        xi = data{i};
        if ~isnumeric(xi)
            error('Cell %d is not numeric.', i);
        end
        cellsOut{i} = xi(:);
        namesOut{i} = sprintf('Var%d', i);
    end
    return

elseif isnumeric(data)
    % numeric matrix: columns are variables, drop NaNs per column
    [nr, nc] = size(data);
    cellsOut = cell(1,nc);
    namesOut = cell(1,nc);
    for i = 1:nc
        cellsOut{i} = data(:, i);
        namesOut{i} = sprintf('Var%d', i);
    end
    return

else
    error('Unsupported input type. Provide table, numeric matrix, or cell array of numeric vectors.');
end
end
