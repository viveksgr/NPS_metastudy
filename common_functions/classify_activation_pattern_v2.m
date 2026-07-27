function results = classify_activation_pattern_v2(pos_act, pos_se, neg_red, neg_se, group, varargin)
% CLASSIFY_ACTIVATION_PATTERN_V2 Analyze a grouped bivariate activation pattern.
%
%   results = classify_activation_pattern_v2(pos_act, pos_se, ...
%       neg_red, neg_se, group)
%
% The function performs only three analyses:
%   1. Bivariate MANOVA on the original Cartesian components [pos_act, neg_red]
%   2. Radial analysis of r = hypot(pos_act, neg_red)
%   3. Angular analysis of theta = atan2(neg_red, pos_act)
%
% INPUTS
%   pos_act : Nx1 estimates for the positive-activation component
%   pos_se  : Nx1 standard errors for pos_act (use zeros or NaNs if unknown)
%   neg_red : Nx1 estimates for the negative-reduction component
%   neg_se  : Nx1 standard errors for neg_red (use zeros or NaNs if unknown)
%   group   : Nx1 numeric, string, cell-string, or categorical group labels
%
% NAME/VALUE OPTIONS
%   'Colors'          : Gx3 RGB matrix (default: lines(G))
%   'MarkerSize'      : scatter marker area (default: 55)
%   'NPermutations'   : radial permutation-test count (default: 10000)
%   'FigureTitle'     : optional overall title
%
% Errors are propagated to polar coordinates with a first-order delta method.
% The inferential tests compare the observed study/subject estimates. The SEs
% and their polar propagation are returned, but are not drawn element-wise.

%% Parse and validate inputs
p = inputParser;
addParameter(p, 'Colors', [], @(x) isempty(x) || ...
    (isnumeric(x) && ismatrix(x) && size(x, 2) == 3));
addParameter(p, 'MarkerSize', 55, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'NPermutations', 10000, @(x) isnumeric(x) && isscalar(x) && x >= 0 && mod(x, 1) == 0);
addParameter(p, 'FigureTitle', '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
opts = p.Results;

pos_act = pos_act(:);
pos_se  = pos_se(:);
neg_red = neg_red(:);
neg_se  = neg_se(:);
group   = group(:);

N = numel(pos_act);
assert(numel(pos_se) == N && numel(neg_red) == N && ...
    numel(neg_se) == N && numel(group) == N, ...
    'All five inputs must have the same number of elements.');
assert(all(pos_se(isfinite(pos_se)) >= 0) && all(neg_se(isfinite(neg_se)) >= 0), ...
    'Standard errors must be nonnegative.');

valid = isfinite(pos_act) & isfinite(neg_red) & ~ismissing(group);
if any(~valid)
    warning('classify_activation_pattern_v2:RowsRemoved', ...
        'Removed %d row(s) with missing estimates or group labels.', sum(~valid));
end
pos_act = pos_act(valid);
pos_se  = pos_se(valid);
neg_red = neg_red(valid);
neg_se  = neg_se(valid);
group   = group(valid);
N = numel(pos_act);

if N < 2
    error('At least two complete observations are required.');
end

if isnumeric(group) || islogical(group)
    group = arrayfun(@(x) sprintf('Group_%g', x), group, 'UniformOutput', false);
end
group_cat   = categorical(group);
group_names = categories(group_cat);
group_idx   = double(group_cat);
n_groups    = numel(group_names);

if n_groups < 2
    error('At least two groups are required for group comparisons.');
end
if any(accumarray(group_idx, 1, [n_groups 1]) < 2)
    warning('classify_activation_pattern_v2:SmallGroup', ...
        'At least one group has fewer than two observations; some tests may be unavailable.');
end

if isempty(opts.Colors)
    colors = lines(n_groups);
else
    if size(opts.Colors, 1) < n_groups
        error('Colors must contain at least one row per group.');
    end
    colors = opts.Colors(1:n_groups, :);
end

%% Cartesian-to-polar transformation and uncertainty propagation
theta_rad = atan2(neg_red, pos_act);
theta_deg = rad2deg(theta_rad);
magnitude = hypot(pos_act, neg_red);

% Delta-method SEs, assuming independent x and y errors.
safe_r = max(magnitude, sqrt(eps));
magnitude_se = sqrt((pos_act ./ safe_r .* pos_se).^2 + ...
                    (neg_red ./ safe_r .* neg_se).^2);
theta_se_rad = sqrt((neg_red ./ safe_r.^2 .* pos_se).^2 + ...
                    (pos_act ./ safe_r.^2 .* neg_se).^2);
theta_se_deg = rad2deg(theta_se_rad);

at_origin = magnitude <= sqrt(eps);
magnitude_se(at_origin) = hypot(pos_se(at_origin), neg_se(at_origin));
theta_se_rad(at_origin) = NaN;
theta_se_deg(at_origin) = NaN;

%% 1. Bivariate MANOVA on [x, y]
p_bivariate_manova = NaN;
manova_dimension = NaN;
manova_stats = [];
try
    [manova_dimension, p_vec, manova_stats] = manova1([pos_act, neg_red], group_idx);
    if ~isempty(p_vec)
        p_bivariate_manova = p_vec(1);
    end
catch ME
    warning('classify_activation_pattern_v2:MANOVAUnavailable', ...
        'Bivariate MANOVA failed: %s', ME.message);
end

%% 2. Radial component analysis
p_radial_anova = NaN;
p_radial_kw = NaN;
p_radial_perm = NaN;
radial_anova_stats = [];

try
    [p_radial_anova, ~, radial_anova_stats] = anova1(magnitude, group_cat, 'off');
catch ME
    warning('classify_activation_pattern_v2:RadialANOVAUnavailable', ...
        'Radial ANOVA failed: %s', ME.message);
end
try
    p_radial_kw = kruskalwallis(magnitude, group_cat, 'off');
catch ME
    warning('classify_activation_pattern_v2:RadialKWUnavailable', ...
        'Radial Kruskal-Wallis test failed: %s', ME.message);
end

if opts.NPermutations > 0
    f_obs = compute_f_stat(magnitude, group_idx, n_groups);
    f_null = NaN(opts.NPermutations, 1);
    for i = 1:opts.NPermutations
        f_null(i) = compute_f_stat(magnitude, group_idx(randperm(N)), n_groups);
    end
    p_radial_perm = (1 + sum(f_null >= f_obs)) / (opts.NPermutations + 1);
end

%% 3. Angular component analysis
p_theta_circ = NaN;
p_theta_anova = NaN;
if exist('circ_wwtest', 'file') == 2
    try
        p_theta_circ = circ_wwtest(theta_rad, group_idx);
    catch ME
        warning('classify_activation_pattern_v2:CircularTestUnavailable', ...
            'CircStat circ_wwtest failed: %s', ME.message);
    end
else
    warning('classify_activation_pattern_v2:CircStatMissing', ...
        'CircStat toolbox not found; circular-test p-value is NaN.');
end
try
    p_theta_anova = anova1(theta_deg, group_cat, 'off');
catch ME
    warning('classify_activation_pattern_v2:AngularANOVAUnavailable', ...
        'Linear angle ANOVA failed: %s', ME.message);
end

%% Group summaries
group_stats = repmat(struct( ...
    'name', '', 'n', 0, ...
    'mean_pos_act', NaN, 'mean_neg_red', NaN, ...
    'mean_magnitude', NaN, 'se_magnitude', NaN, ...
    'circular_mean_theta_deg', NaN, 'resultant_length', NaN), n_groups, 1);

for g = 1:n_groups
    idx = group_idx == g;
    zg = mean(exp(1i * theta_rad(idx)));
    group_stats(g).name = group_names{g};
    group_stats(g).n = sum(idx);
    group_stats(g).mean_pos_act = mean(pos_act(idx));
    group_stats(g).mean_neg_red = mean(neg_red(idx));
    group_stats(g).mean_magnitude = mean(magnitude(idx));
    group_stats(g).se_magnitude = std(magnitude(idx)) / sqrt(sum(idx));
    group_stats(g).circular_mean_theta_deg = rad2deg(angle(zg));
    group_stats(g).resultant_length = abs(zg);
end

%% Plot the three requested analyses
fig = figure('Color', 'w', 'Position', [60 180 1500 470]);
t = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: original bivariate components
ax_xy = nexttile(t, 1);
hold(ax_xy, 'on');
box(ax_xy, 'on');
xline(ax_xy, 0, '-', 'Color', [0.65 0.65 0.65], 'HandleVisibility', 'off');
yline(ax_xy, 0, '-', 'Color', [0.65 0.65 0.65], 'HandleVisibility', 'off');
for g = 1:n_groups
    idx = find(group_idx == g);
    scatter(ax_xy, pos_act(idx), neg_red(idx), opts.MarkerSize, colors(g, :), ...
        'filled', 'MarkerFaceAlpha', 0.85, 'DisplayName', group_names{g});
end

% If MANOVA is significant, show the first canonical separation axis
% through the overall centroid (as in the original function).
if p_bivariate_manova < 0.05 && manova_dimension >= 1 && ...
        ~isempty(manova_stats) && isfield(manova_stats, 'eigenvec')
    v = manova_stats.eigenvec(:, 1);
    if numel(v) >= 2 && all(isfinite(v(1:2))) && norm(v(1:2)) > 0
        v = v(1:2) / norm(v(1:2));
        centroid = [mean(pos_act), mean(neg_red)];
        xy_limits = axis(ax_xy);
        line_radius = hypot(diff(xy_limits(1:2)), diff(xy_limits(3:4)));
        plot(ax_xy, centroid(1) + [-line_radius line_radius] * v(1), ...
            centroid(2) + [-line_radius line_radius] * v(2), ...
            'k:', 'LineWidth', 1.6, 'HandleVisibility', 'off');
        axis(ax_xy, xy_limits);
    end
end
xlabel(ax_xy, 'Positive activation');
ylabel(ax_xy, 'Negative reduction');
title(ax_xy, sprintf('Bivariate MANOVA, p = %s', format_p(p_bivariate_manova)));
axis(ax_xy, 'square');
legend(ax_xy, 'Location', 'best', 'Interpreter', 'none');

% Panel 2: radial magnitude by group
ax_r = nexttile(t, 2);
hold(ax_r, 'on');
box(ax_r, 'on');
boxplot(ax_r, magnitude, group_cat, 'Colors', colors, 'Symbol', '', ...
    'GroupOrder', group_names);
jitter = deterministic_jitter(group_idx);
for g = 1:n_groups
    idx = group_idx == g;
    scatter(ax_r, group_idx(idx) + jitter(idx), magnitude(idx), ...
        opts.MarkerSize * 0.65, colors(g, :), 'filled', ...
        'MarkerFaceAlpha', 0.7, 'HandleVisibility', 'off');
    scatter(ax_r, g, group_stats(g).mean_magnitude, opts.MarkerSize * 1.5, ...
        colors(g, :), 'd', 'filled', 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');
end
ylabel(ax_r, 'Radius, r');
title(ax_r, sprintf('Radial: ANOVA %s, KW %s, perm %s', ...
    format_p(p_radial_anova), format_p(p_radial_kw), format_p(p_radial_perm)));
xticklabels(ax_r, group_names);
xtickangle(ax_r, 25);

% Panel 3: theta-only polar plot. Radius is fixed to one so that radial
% magnitude is represented only in Panel 2.
polar_parent = nexttile(t, 3);
polar_position = polar_parent.Position;
delete(polar_parent);
ax_polar = polaraxes('Parent', fig, 'Position', polar_position);
hold(ax_polar, 'on');
ax_polar.ThetaZeroLocation = 'right';
ax_polar.ThetaDir = 'counterclockwise';

for g = 1:n_groups
    idx = group_idx == g;
    polarscatter(ax_polar, theta_rad(idx), ones(sum(group_idx == g), 1), opts.MarkerSize, ...
        colors(g, :), 'filled', 'MarkerFaceAlpha', 0.30, ...
        'MarkerEdgeAlpha', 0.30, ...
        'DisplayName', group_names{g});
    polarscatter(ax_polar, deg2rad(group_stats(g).circular_mean_theta_deg), 1, ...
        opts.MarkerSize * 1.7, colors(g, :), 'd', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end
ax_polar.RLim = [0 1.12];
ax_polar.RTick = 1;
ax_polar.RTickLabel = {''};
title(ax_polar, sprintf('Angular: circular %s, linear ANOVA %s', ...
    format_p(p_theta_circ), format_p(p_theta_anova)));
legend(ax_polar, 'Location', 'bestoutside', 'Interpreter', 'none');

if strlength(string(opts.FigureTitle)) > 0
    title(t, char(opts.FigureTitle), 'Interpreter', 'none');
end

%% Console summary
fprintf('\n=== Activation Pattern Analysis v2 ===\n');
fprintf('Bivariate MANOVA [x,y] : p = %s\n', format_p(p_bivariate_manova));
fprintf('Radial magnitude       : ANOVA p = %s | KW p = %s | permutation p = %s\n', ...
    format_p(p_radial_anova), format_p(p_radial_kw), format_p(p_radial_perm));
fprintf('Angular direction      : circular p = %s | linear ANOVA p = %s\n\n', ...
    format_p(p_theta_circ), format_p(p_theta_anova));

%% Pack results
results.pos_act = pos_act;
results.pos_se = pos_se;
results.neg_red = neg_red;
results.neg_se = neg_se;
results.group = group_cat;
results.group_names = group_names;
results.colors = colors;
results.magnitude = magnitude;
results.magnitude_se = magnitude_se;
results.theta_rad = theta_rad;
results.theta_deg = theta_deg;
results.theta_se_rad = theta_se_rad;
results.theta_se_deg = theta_se_deg;
results.p_bivariate_manova = p_bivariate_manova;
results.manova_dimension = manova_dimension;
results.manova_stats = manova_stats;
results.p_radial_anova = p_radial_anova;
results.p_radial_kw = p_radial_kw;
results.p_radial_perm = p_radial_perm;
results.radial_anova_stats = radial_anova_stats;
results.p_theta_circ = p_theta_circ;
results.p_theta_anova = p_theta_anova;
results.group_stats = group_stats;
results.figure = fig;
results.axes = struct('bivariate', ax_xy, 'radial', ax_r, 'angular_polar', ax_polar);

end

function f = compute_f_stat(y, group_idx, n_groups)
grand_mean = mean(y);
ss_between = 0;
ss_within = 0;
for g = 1:n_groups
    yg = y(group_idx == g);
    if isempty(yg)
        f = NaN;
        return;
    end
    ss_between = ss_between + numel(yg) * (mean(yg) - grand_mean)^2;
    ss_within = ss_within + sum((yg - mean(yg)).^2);
end
df_between = n_groups - 1;
df_within = numel(y) - n_groups;
if df_between <= 0 || df_within <= 0
    f = NaN;
else
    f = (ss_between / df_between) / (ss_within / df_within + eps);
end
end

function jitter = deterministic_jitter(group_idx)
jitter = zeros(size(group_idx));
groups = unique(group_idx(:))';
for g = groups
    idx = find(group_idx == g);
    if numel(idx) > 1
        jitter(idx) = linspace(-0.10, 0.10, numel(idx));
    end
end
end

function txt = format_p(p)
if isempty(p) || ~isfinite(p)
    txt = 'n/a';
elseif p < 0.001
    txt = '<.001';
else
    txt = sprintf('%.3f', p);
end
end
