function [fig, order, S] = spectralReorderAndPlot(M, varargin)
% SPECTRALREORDERANDPLOT  Spectral reorder columns and plot heatmap
%
%   [fig, order, S] = spectralReorderAndPlot(M)
%   [fig, order, S] = spectralReorderAndPlot(M, 'Labels', labels, 'Affinity', 'corr', ...)
%
% INPUT
%   M         - nROI x nInterventions numeric matrix (can contain NaNs)
%
% NAME-VALUE OPTIONS
%   'Labels'      - cellstr of length nInterventions for x tick labels (default {})
%   'Affinity'    - 'corr' (default) | 'abs_corr' | 'euclid_gauss'
%   'GaussianSigma' - scalar sigma for 'euclid_gauss' (default = median pairwise dist)
%   'Clim'        - 2-element [min max] color limits (default symmetric about 0)
%   'Title'       - figure title string (default '')
%   'Colormap'    - colormap handle or matrix (default blue-white-red)
%   'FigureSize'  - [w h] in pixels (default [900 500])
%
% OUTPUT
%   fig   - figure handle
%   order - permutation of columns (1..nInterventions) after reordering
%   S     - column similarity (nInterventions x nInterventions)
%
% EXAMPLE
%   % M is nROI x nInt
%   [fig, ord] = spectralReorderAndPlot(M, 'Labels', interventionNames);
%

% ---------------- parse inputs ----------------
p = inputParser;
addParameter(p,'Labels',{}, @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(p,'Affinity','corr', @(s) ismember(s,{'corr','abs_corr','euclid_gauss'}));
addParameter(p,'GaussianSigma',[], @isnumeric);
addParameter(p,'Clim',[], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(p,'Title','', @ischar);
addParameter(p,'Colormap', [], @(x) isempty(x) || (isnumeric(x) || isa(x,'function_handle')));
addParameter(p,'FigureSize',[900 500], @(x)isnumeric(x) && numel(x)==2);
parse(p, varargin{:});
labels = p.Results.Labels;
affinity_mode = p.Results.Affinity;
sigma = p.Results.GaussianSigma;
clim_in = p.Results.Clim;
figTitle = p.Results.Title;
cmap_in = p.Results.Colormap;
figSize = p.Results.FigureSize;

[nROI, nInt] = size(M);
if nInt < 2, error('Need at least 2 interventions (columns).'); end

% ---------------- compute column similarity matrix S ----------------
switch affinity_mode
    case 'corr'
        S = corr(M, 'rows', 'pairwise');   % nInt x nInt (columns correlations)
        S(isnan(S)) = 0;
        S(S < 0) = 0;                       % clip negatives -> nonnegative affinity
    case 'abs_corr'
        S = abs(corr(M, 'rows', 'pairwise'));
        S(isnan(S)) = 0;
    case 'euclid_gauss'
        % compute pairwise euclidean distances between columns (ignoring NaNs pairwise)
        D = zeros(nInt);
        for i=1:nInt
            for j=i+1:nInt
                xi = M(:,i); xj = M(:,j);
                valid = ~isnan(xi) & ~isnan(xj);
                if sum(valid) < 1
                    dij = inf;
                else
                    dij = norm(xi(valid)-xj(valid));
                end
                D(i,j) = dij; D(j,i) = dij;
            end
        end
        if isempty(sigma)
            finiteD = D(isfinite(D) & D>0);
            if isempty(finiteD), sigma = 1; else sigma = median(finiteD(:)); end
        end
        S = exp(-(D.^2) / (2 * (sigma^2 + eps)));
        S(~isfinite(S)) = 0;
    otherwise
        error('Unknown affinity mode.');
end

% ensure diagonal is 0 (not used)
S(1:nInt+1:end) = 0;

% ---------------- spectral reordering (Fiedler) ----------------
order = 1:nInt;
try
    % degree and normalized Laplacian
    d = sum(S,2);
    DinvSqrt = diag(1./sqrt(d + eps));
    Lsym = eye(nInt) - DinvSqrt * S * DinvSqrt;   % normalized Laplacian
    % eigen decomposition (symmetric)
    opts.issym = true; opts.isreal = true;
    [V, E] = eig((Lsym + Lsym')/2);
    eigvals = diag(E);
    [eigvals_sorted, idx] = sort(eigvals, 'ascend');
    % find second smallest eigenvector (Fiedler). First eigenvector ~ constant.
    if numel(eigvals_sorted) < 2
        error('Not enough eigenvalues.');
    end
    fidx = idx(2);                % index into V for Fiedler
    fvec = V(:, fidx);
    % sort by Fiedler component
    [~, order] = sort(fvec, 'ascend');
    % if it leads to trivial order (all equal), fallback
    if numel(unique(fvec)) == 1
        error('Fiedler vector constant -> fallback');
    end
catch
    % fallback: hierarchical clustering ordering
    try
        Y = pdist(S, 'euclidean');
        Z = linkage(Y, 'average');
        order = optimalleaforder(Z, Y);
    catch
        % last fallback: simple correlation-based sorting by first PC of S
        [~,~,Vpc] = svd(S, 'econ');
        pc1 = Vpc(:,1);
        [~, order] = sort(pc1, 'descend');
    end
end

% ---------------- plotting ----------------
fig = figure('Color','w', 'Position', [100 100 figSize]);
% imagesc expects rows = y, cols = x; we want y=ROI, x=intervention:
% display transposed so each row is an intervention: imagesc(1:nInt, 1:nROI, M(:,order)')
if isempty(clim_in)
    mx = max(abs(M(:)));
    clim = [-mx mx];
else
    clim = clim_in;
end

% default diverging cmap if none provided
if isempty(cmap_in)
    ncol = 256;
    half = floor(ncol/2);
    c1 = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];       % blue->white
    c2 = [ones(ncol-half,1), linspace(1,0,ncol-half)', linspace(1,0,ncol-half)']; % white->red
    cmap = [c1; c2];
else
    if isa(cmap_in, 'function_handle'), cmap = cmap_in(256); else cmap = cmap_in; end
end

h = imagesc(1:nInt, 1:nROI, M(:,order)', clim);
colormap(cmap); colorbar;
set(gca, 'YDir','normal');
xlabel('Intervention (reordered)');
ylabel('ROI');
title(figTitle, 'Interpreter','none');

% xtick labels if provided
if ~isempty(labels)
    labels = cellstr(labels);
    if numel(labels) == nInt
        set(gca, 'XTick', 1:nInt, 'XTickLabel', labels(order), 'XTickLabelRotation', 45, 'TickLabelInterpreter','none');
    else
        warning('Labels length mismatch; skipping xticklabels.');
    end
else
    set(gca, 'XTick', 1:nInt);
end

% improve visual (tight, grid)
axis tight; box on;

% ensure NaNs appear distinct (map NaN to gray)
% Create an alpha mask so NaNs are transparent and set axes background to light gray
nanMask = isnan(M(:,order)');
if any(nanMask(:))
    set(h, 'AlphaData', ~nanMask);
    set(gca, 'Color', [0.9 0.9 0.9]); % background gray for NaNs
end

end
