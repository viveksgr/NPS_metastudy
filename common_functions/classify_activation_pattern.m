function results = classify_activation_pattern(pos_act, neg_red, group)
% Classifies subjects by response pattern and tests whether group membership
% predicts (a) radial strength (magnitude), (b) angular signature (theta),
% (c) joint r-theta signature (MANOVA), (d) axial angular signature (folded),
% and (e) quadrant proportions (chi-square).
%
% Sign convention: both pos_act and neg_red are expected to be POSITIVE when
% the intervention "works" (pos ROIs activate, neg ROIs reduce - multiplied by
% -1 upstream). Data can span all four quadrants.
%
% Zones (by theta = atan2d(neg_red, pos_act), in degrees) -- full quadrants:
%   both             : theta in [   0,  90)   -> Q1
%   suppression_only : theta in [  90, 180]   -> Q2
%   neither          : theta in [-180, -90)   -> Q3
%   activation_only  : theta in [ -90,   0)   -> Q4
%
% Inputs:
%   pos_act : Nx1 mean activation in positive ROIs (positive = good)
%   neg_red : Nx1 mean reduction in negative ROIs  (positive = good)
%   group   : Nx1 group labels (cell array of strings OR numeric vector)
%
% Outputs (results struct):
%   .theta_deg, .theta_rad, .magnitude, .balance, .type
%   .p_radial_anova, .p_radial_kw, .p_radial_perm    (magnitude tests)
%   .p_theta_circ, .p_theta_anova                    (angular tests, unfolded)
%   .p_theta_circ_axial, .p_theta_anova_axial        (angular tests, folded mod pi)
%   .p_rtheta_manova                                 (MANOVA on [r, theta])
%   .p_chi2                                          (chi-square on quadrant %)
%   .posthoc_radial                                  (Tukey HSD on magnitude)
%   .pct_matrix, .group_stats, .group_names

%% --- Parse group input (canonical order via categorical) ---
N = numel(pos_act);
assert(numel(neg_red) == N && numel(group) == N, 'All inputs must have the same length.');

if isnumeric(group)
    group = arrayfun(@(x) sprintf('Group_%g', x), group, 'UniformOutput', false);
end
group_cat   = categorical(group);
group_names = categories(group_cat);          % canonical (sorted) order
group_idx   = double(group_cat);              % 1..n_groups matching categories()
n_groups    = numel(group_names);
colors      = lines(n_groups);                % color g -> group_names{g}

%% --- Core geometry ---
theta_deg = atan2d(neg_red, pos_act);          % [-180, 180]
theta_rad = atan2(neg_red, pos_act);
magnitude = sqrt(pos_act.^2 + neg_red.^2);

% Axial-folded versions (treat theta and theta-pi as equivalent)
theta_deg_axial = mod(theta_deg + 90, 180) - 90;   % [-90, 90]
theta_rad_axial = 2 * theta_rad;                   % doubled-angle convention

% Balance (Q1 only)
in_q1 = pos_act > 0 & neg_red > 0;
balance = zeros(N, 1);
balance(in_q1) = 2 .* min(pos_act(in_q1), neg_red(in_q1)) ./ ...
                 (pos_act(in_q1) + neg_red(in_q1) + eps);

%% --- Categorical classification (full quadrants) ---
type = repmat({'neither'}, N, 1);
type(theta_deg >= -90  & theta_deg <   0)  = {'activation_only'};   % Q4
type(theta_deg >=   0  & theta_deg <  90)  = {'both'};              % Q1
type(theta_deg >=  90  & theta_deg <= 180) = {'suppression_only'};  % Q2
% theta_deg in [-180, -90) remains 'neither'                         % Q3

%% --- Radial (magnitude) tests ---
[p_radial_anova, ~, stats_radial] = anova1(magnitude, group_cat, 'off');
[p_radial_kw,    ~]               = kruskalwallis(magnitude, group_cat, 'off');

f_obs_r  = compute_f_stat(magnitude, group_idx, n_groups);
n_perm   = 10000;
f_null_r = zeros(n_perm, 1);
for i = 1:n_perm
    f_null_r(i) = compute_f_stat(magnitude, group_idx(randperm(N)), n_groups);
end
p_radial_perm = mean(f_null_r >= f_obs_r);

c_r = multcompare(stats_radial, 'Display', 'off');
posthoc_radial = table(group_names(c_r(:,1)), group_names(c_r(:,2)), c_r(:,6), ...
    'VariableNames', {'Group1', 'Group2', 'p_tukey'});

%% --- Angular (theta) tests: unfolded and axial (folded mod pi) ---
has_circstat = exist('circ_wwtest', 'file') == 2;
p_theta_circ       = NaN;
p_theta_circ_axial = NaN;
if has_circstat
    try, p_theta_circ       = circ_wwtest(theta_rad,       group_idx); catch, end
    try, p_theta_circ_axial = circ_wwtest(theta_rad_axial, group_idx); catch, end
else
    warning('CircStat toolbox not found; circ_wwtest p-values set to NaN.');
end
[p_theta_anova,       ~] = anova1(theta_deg,       group_cat, 'off');
[p_theta_anova_axial, ~] = anova1(theta_deg_axial, group_cat, 'off');

%% --- MANOVA on [r, theta] ---
try
    [~, p_manova_vec] = manova1([magnitude, theta_deg], group_idx);
    if isempty(p_manova_vec)
        p_rtheta_manova = NaN;
    else
        p_rtheta_manova = p_manova_vec(1);
    end
catch ME
    warning('manova1 on [r,theta] failed: %s', ME.message);
    p_rtheta_manova = NaN;
end

%% --- MANOVA on [pos_act, neg_red] (for subplot 1,1) ---
try
    [d_xy, p_xy_vec, stats_xy] = manova1([pos_act, neg_red], group_idx);
    if isempty(p_xy_vec)
        p_manova_xy = NaN;
    else
        p_manova_xy = p_xy_vec(1);
    end
catch ME
    warning('manova1 on [pos,neg] failed: %s', ME.message);
    d_xy = 0; p_manova_xy = NaN; stats_xy = [];
end

%% --- Chi-square on quadrant proportions ---
type_names_order = {'activation_only', 'both', 'suppression_only', 'neither'};
type_idx = zeros(N, 1);
for k = 1:numel(type_names_order)
    type_idx(strcmp(type, type_names_order{k})) = k;
end
try
    [~, ~, p_chi2] = crosstab(group_idx, type_idx);
catch
    p_chi2 = NaN;
end

%% --- Per-group summary ---
group_stats = struct();
for g = 1:n_groups
    idx = group_idx == g;
    group_stats(g).name           = group_names{g};
    group_stats(g).n              = sum(idx);
    group_stats(g).mean_magnitude = mean(magnitude(idx));
    group_stats(g).mean_theta     = mean(theta_deg(idx));
    group_stats(g).mean_balance   = mean(balance(idx));
    group_stats(g).pct_act_only   = 100 * mean(strcmp(type(idx), 'activation_only'));
    group_stats(g).pct_both       = 100 * mean(strcmp(type(idx), 'both'));
    group_stats(g).pct_sup_only   = 100 * mean(strcmp(type(idx), 'suppression_only'));
    group_stats(g).pct_neither    = 100 * mean(strcmp(type(idx), 'neither'));
end

pct_matrix = zeros(n_groups, 4);
for g = 1:n_groups
    pct_matrix(g, :) = [group_stats(g).pct_act_only, group_stats(g).pct_both, ...
                        group_stats(g).pct_sup_only, group_stats(g).pct_neither];
end

%% --- Plotting: 2x3 layout ---
figure('Position', [60 80 1600 950], 'Color', 'w');

% ---- Panel 1: 2D scatter (x,y) with quadrant zones ----
subplot(2, 3, 1); hold on; box on;
xmin = min([0; pos_act]) * 1.15;  xmax = max([0; pos_act]) * 1.15;
ymin = min([0; neg_red]) * 1.15;  ymax = max([0; neg_red]) * 1.15;
R = 1.5 * max(abs([xmin xmax ymin ymax]));

draw_wedge(   0,  90, R, [0.75 1.0 0.75]);  % both            (Q1)
draw_wedge(  90, 180, R, [0.75 0.75 1.0]);  % suppression_only(Q2)
draw_wedge(-180, -90, R, [0.80 0.80 0.80]); % neither         (Q3)
draw_wedge( -90,   0, R, [1.0 0.75 0.75]);  % activation_only (Q4)

plot([xmin xmax], [0 0], 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot([0 0], [ymin ymax], 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');

for g = 1:n_groups
    idx = group_idx == g;
    scatter(pos_act(idx), neg_red(idx), 60, colors(g,:), 'filled', ...
            'DisplayName', group_names{g}, 'MarkerFaceAlpha', 0.85);
end

% Plot canonical separation axes through centroid if MANOVA indicates separation
if d_xy >= 1 && ~isempty(stats_xy) && isfield(stats_xy, 'eigenvec')
    centroid = [mean(pos_act), mean(neg_red)];
    for k = 1:d_xy
        v = stats_xy.eigenvec(:, k);
        if norm(v) > 0
            v = v / norm(v);
            plot(centroid(1) + [-R R]*v(1), centroid(2) + [-R R]*v(2), ...
                 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
end

xlabel('Positive activation'); ylabel('Negative reduction (sign-flipped)');
if isnan(p_manova_xy)
    title('Activation pattern (x, y)');
else
    title(sprintf('Activation pattern (x, y)   MANOVA d=%d  p=%.3f', d_xy, p_manova_xy));
end
legend('Location', 'best', 'FontSize', 8);
xlim([xmin xmax]); ylim([ymin ymax]); axis square;

% ---- Panel 2: Radial (magnitude) by group ----
subplot(2, 3, 2); hold on; box on;
boxplot(magnitude, group_cat, 'Colors', colors, 'Symbol', '', 'GroupOrder', group_names);
for i = 1:N
    xpos = group_idx(i) + 0.12*(rand - 0.5);
    scatter(xpos, magnitude(i), 35, colors(group_idx(i),:), 'filled', ...
            'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
end
ylabel('Magnitude  \surd(pos^2 + neg^2)');
title(sprintf('Radial strength   ANOVA p=%.3f   KW p=%.3f   Perm p=%.3f', ...
      p_radial_anova, p_radial_kw, p_radial_perm));
xticklabels(group_names); xtickangle(30);

% ---- Panel 3: Angular (theta) by group, unfolded [-180, 180] ----
subplot(2, 3, 3); hold on; box on;
xl = [0.5, n_groups + 0.5];
draw_theta_bands(xl, false);
for yv = [-90, 0, 90]
    yline(yv, '-', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off');
end
boxplot(theta_deg, group_cat, 'Colors', colors, 'Symbol', '', 'GroupOrder', group_names);
for i = 1:N
    xpos = group_idx(i) + 0.12*(rand - 0.5);
    scatter(xpos, theta_deg(i), 35, colors(group_idx(i),:), 'filled', ...
            'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
end
ylabel('\theta (deg)');
if isnan(p_theta_circ)
    title(sprintf('Angular signature   ANOVA p=%.3f   (circ\\_ww n/a)', p_theta_anova));
else
    title(sprintf('Angular signature   circ\\_ww p=%.3f   ANOVA p=%.3f', ...
          p_theta_circ, p_theta_anova));
end
xticklabels(group_names); xtickangle(30);
xlim(xl); ylim([-180 180]); yticks(-180:45:180);

% ---- Panel 4: r vs theta scatter + MANOVA1 ----
subplot(2, 3, 4); hold on; box on;

% Shaded vertical bands mirroring quadrant coloring along theta axis
yl4 = [0, max(magnitude) * 1.15];
draw_theta_vbands(yl4);
for xv = [-90, 0, 90]
    xline(xv, '-', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off');
end

for g = 1:n_groups
    idx = group_idx == g;
    scatter(theta_deg(idx), magnitude(idx), 60, colors(g,:), 'filled', ...
            'DisplayName', group_names{g}, 'MarkerFaceAlpha', 0.85);
end
xlabel('\theta (deg)'); ylabel('r = magnitude');
if isnan(p_rtheta_manova)
    title('(r, \theta)   MANOVA unavailable');
else
    title(sprintf('(r, \\theta)   MANOVA p=%.3f', p_rtheta_manova));
end
xlim([-180 180]); xticks(-180:45:180); ylim(yl4);
legend('Location', 'best', 'FontSize', 8);

% ---- Panel 5: Axial (folded) theta by group ----
subplot(2, 3, 5); hold on; box on;
xl = [0.5, n_groups + 0.5];
draw_theta_bands(xl, true);   % two folded bands
yline(0, '-', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off');

boxplot(theta_deg_axial, group_cat, 'Colors', colors, 'Symbol', '', 'GroupOrder', group_names);
for i = 1:N
    xpos = group_idx(i) + 0.12*(rand - 0.5);
    scatter(xpos, theta_deg_axial(i), 35, colors(group_idx(i),:), 'filled', ...
            'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
end

ylabel('\theta (deg, folded mod \pi)');
if isnan(p_theta_circ_axial)
    title(sprintf('Axial (folded)   ANOVA p=%.3f   (circ\\_ww n/a)', p_theta_anova_axial));
else
    title(sprintf('Axial (folded)   circ\\_ww p=%.3f   ANOVA p=%.3f', ...
          p_theta_circ_axial, p_theta_anova_axial));
end
xticklabels(group_names); xtickangle(30);
xlim(xl); ylim([-90 90]); yticks(-90:30:90);

% ---- Panel 6: Stacked bar of quadrant proportions ----
subplot(2, 3, 6);
b = bar(pct_matrix, 'stacked');
b(1).FaceColor = [0.90 0.50 0.50];  % activation_only
b(2).FaceColor = [0.50 0.85 0.50];  % both
b(3).FaceColor = [0.50 0.50 0.90];  % suppression_only
b(4).FaceColor = [0.70 0.70 0.70];  % neither
xticklabels(group_names); xtickangle(30);
ylabel('% of studies');
if isnan(p_chi2)
    title('Quadrant proportions');
else
    title(sprintf('Quadrant proportions   \\chi^2 p=%.3f', p_chi2));
end
legend({'Act only','Both','Sup only','Neither'}, 'Location','eastoutside','FontSize',8);
ylim([0 105]); box on;

%% --- Console summary ---
fprintf('\n=== Activation Pattern Analysis ===\n');
fprintf('Radial (magnitude)  — ANOVA p=%.4f | KW p=%.4f | Perm p=%.4f\n', ...
        p_radial_anova, p_radial_kw, p_radial_perm);
if isnan(p_theta_circ)
    fprintf('Angular (unfolded)  — circ_ww n/a | ANOVA p=%.4f\n', p_theta_anova);
else
    fprintf('Angular (unfolded)  — circ_ww p=%.4f | ANOVA p=%.4f\n', p_theta_circ, p_theta_anova);
end
if isnan(p_theta_circ_axial)
    fprintf('Angular (axial)     — circ_ww n/a | ANOVA p=%.4f\n', p_theta_anova_axial);
else
    fprintf('Angular (axial)     — circ_ww p=%.4f | ANOVA p=%.4f\n', p_theta_circ_axial, p_theta_anova_axial);
end
fprintf('(r, theta)          — MANOVA p=%.4f\n', p_rtheta_manova);
fprintf('Quadrant proportions— chi2 p=%.4f\n\n', p_chi2);

fprintf('%-20s %4s %8s %8s %8s %8s %8s %8s\n', ...
        'Group','N','|v|','theta','Act%','Both%','Sup%','None%');
fprintf('%s\n', repmat('-', 1, 80));
for g = 1:n_groups
    fprintf('%-20s %4d %8.3f %8.1f %7.1f%% %7.1f%% %7.1f%% %7.1f%%\n', ...
        group_stats(g).name, group_stats(g).n, ...
        group_stats(g).mean_magnitude, group_stats(g).mean_theta, ...
        group_stats(g).pct_act_only, group_stats(g).pct_both, ...
        group_stats(g).pct_sup_only, group_stats(g).pct_neither);
end

fprintf('\nPost-hoc pairwise on magnitude (Tukey HSD):\n');
for i = 1:height(posthoc_radial)
    sig = '';
    if posthoc_radial.p_tukey(i) < 0.001, sig = '***';
    elseif posthoc_radial.p_tukey(i) < 0.01, sig = '**';
    elseif posthoc_radial.p_tukey(i) < 0.05, sig = '*';
    end
    fprintf('  %-18s vs %-18s  p=%.4f %s\n', ...
        posthoc_radial.Group1{i}, posthoc_radial.Group2{i}, posthoc_radial.p_tukey(i), sig);
end

%% --- Pack results ---
results.theta_deg           = theta_deg;
results.theta_rad           = theta_rad;
results.theta_deg_axial     = theta_deg_axial;
results.magnitude           = magnitude;
results.balance             = balance;
results.type                = type;
results.p_radial_anova      = p_radial_anova;
results.p_radial_kw         = p_radial_kw;
results.p_radial_perm       = p_radial_perm;
results.p_theta_circ        = p_theta_circ;
results.p_theta_anova       = p_theta_anova;
results.p_theta_circ_axial  = p_theta_circ_axial;
results.p_theta_anova_axial = p_theta_anova_axial;
results.p_rtheta_manova     = p_rtheta_manova;
results.p_manova_xy         = p_manova_xy;
results.d_manova_xy         = d_xy;
results.p_chi2              = p_chi2;
results.posthoc_radial      = posthoc_radial;
results.pct_matrix          = pct_matrix;
results.group_stats         = group_stats;
results.group_names         = group_names;

end

%% --- Helper: angular wedge for scatter panel ---
function draw_wedge(a1_deg, a2_deg, R, color)
n = 40;
t = deg2rad(linspace(a1_deg, a2_deg, n));
x = [0, R*cos(t), 0];
y = [0, R*sin(t), 0];
fill(x, y, color, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

%% --- Helper: horizontal theta bands behind boxplots (panels 3, 5) ---
function draw_theta_bands(xl, folded)
if ~folded
    patch([xl fliplr(xl)], [-180 -180 -90 -90], [0.80 0.80 0.80], ...
          'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
    patch([xl fliplr(xl)], [ -90  -90   0   0], [1.00 0.75 0.75], ...
          'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
    patch([xl fliplr(xl)], [   0    0  90  90], [0.75 1.00 0.75], ...
          'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
    patch([xl fliplr(xl)], [  90   90 180 180], [0.75 0.75 1.00], ...
          'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
else
    patch([xl fliplr(xl)], [-90 -90  0  0], [0.75 0.60 0.85], ...
          'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
    patch([xl fliplr(xl)], [  0   0 90 90], [0.60 0.85 0.75], ...
          'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
end
end

%% --- Helper: vertical theta bands behind r-theta scatter (panel 4) ---
function draw_theta_vbands(yl)
% Rectangles along theta axis for the 4 quadrants
patch([-180 -90 -90 -180], [yl(1) yl(1) yl(2) yl(2)], [0.80 0.80 0.80], ...
      'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
patch([ -90   0   0  -90], [yl(1) yl(1) yl(2) yl(2)], [1.00 0.75 0.75], ...
      'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
patch([   0  90  90    0], [yl(1) yl(1) yl(2) yl(2)], [0.75 1.00 0.75], ...
      'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
patch([  90 180 180   90], [yl(1) yl(1) yl(2) yl(2)], [0.75 0.75 1.00], ...
      'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
end

%% --- Helper: F-statistic for permutation test ---
function f = compute_f_stat(y, group_idx, n_groups)
grand_mean = mean(y);
N = numel(y);
ss_between = 0;
ss_within  = 0;
for g = 1:n_groups
    yg = y(group_idx == g);
    ng = numel(yg);
    ss_between = ss_between + ng * (mean(yg) - grand_mean)^2;
    ss_within  = ss_within  + sum((yg - mean(yg)).^2);
end
df_between = n_groups - 1;
df_within  = N - n_groups;
f = (ss_between / df_between) / (ss_within / df_within + eps);
end
