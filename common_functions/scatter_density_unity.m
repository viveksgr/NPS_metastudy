function [fig, stats] = scatter_density_unity(x, y, varargin)
% SCATTER_DENSITY_UNITY  Bivariate histogram with unity line and test.
%
%   [fig, stats] = scatter_density_unity(x, y)
%   [fig, stats] = scatter_density_unity(x, y, 'xlabel', 'A', 'ylabel', 'B', 'nbins', 200)
%
% Plots a 2-D histogram (imagesc) of x vs y, overlays the unity line,
% and tests whether the distribution sits above or below unity
% (sign-rank test on y - x). The p-value and direction are shown in the title.
%
% INPUTS
%   x, y      - numeric vectors of equal length
%
% OPTIONAL NAME/VALUE
%   'xlabel'  - x-axis label (default 'X')
%   'ylabel'  - y-axis label (default 'Y')
%   'nbins'   - number of bins per axis (default 200)
%
% OUTPUTS
%   fig       - figure handle
%   stats     - struct with fields: p, signedrank, median_diff, direction

p = inputParser;
addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addParameter(p, 'xlabel', 'X', @ischar);
addParameter(p, 'ylabel', 'Y', @ischar);
addParameter(p, 'nbins', 200, @isscalar);
parse(p, x, y, varargin{:});

xlab = p.Results.xlabel;
ylab = p.Results.ylabel;
nbins = p.Results.nbins;

x = x(:);
y = y(:);
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

% Bivariate histogram (shared axis range, integer-aligned edges)
lo = floor(min([x; y]));
hi = ceil(max([x; y]));
step = (hi - lo) / nbins;
xedges = lo:step:hi;
yedges = xedges;
C = histcounts2(x, y, xedges, yedges);

% Centers for axis labels
xc = xedges(1:end-1) + diff(xedges)/2;
yc = yedges(1:end-1) + diff(yedges)/2;

fig = figure('Position', [200 200 560 480]);
imagesc(xc, yc, C');
axis xy equal tight;
set(gca, 'XGrid', 'off', 'YGrid', 'off', 'GridLineStyle', 'none', ...
         'MinorGridLineStyle', 'none', 'Layer', 'top', 'TickDir', 'out');
colormap(turbo);
colorbar;
hold on;

lims = [min([x; y]), max([x; y])];
plot(lims, lims, 'w--', 'LineWidth', 1.5);
hold off;

xlabel(xlab);
ylabel(ylab);

% Wilcoxon signed-rank test on y - x
d = y - x;
[pval, ~, sr] = signrank(d);
med_diff = median(d);
if med_diff > 0
    dir_str = 'above';
else
    dir_str = 'below';
end

title(sprintf('Median diff = %.3g (%s unity), p = %.3g', med_diff, dir_str, pval));

stats.p = pval;
stats.signedrank = sr;
stats.median_diff = med_diff;
stats.direction = dir_str;

end
