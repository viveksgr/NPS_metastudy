function out = run_readout(varargin)
% RUN_READOUT  Spec Step A (§3): estimate the neural->behavior readout map w.
%
% w is the WITHIN-STUDY, within-subject slope of relief (=-pain) on the
% neural coordinates, NOT the cross-study ecological correlation. Fit an LMM
% pooled across all subjects with a study random intercept (minimum), and,
% where estimable, random slopes. Also contrasts w against the ecological
% (study-mean) slope, and characterizes the low-rank structure + null space.
%
% This is estimable despite the study==intervention confound: it concerns the
% neural->behavior mapping, not between-group contrasts.
%
% Usage:
%   out = run_readout();                          % COORDS = {C1,C2,C3}
%   out = run_readout('COORDS',{'NPS','SIIPS'});
%   out = run_readout('COORDS',{'NPS','SIIPS','C1','C2','C3'});

p = inputParser;
addParameter(p,'MAT','datamat.mat');
addParameter(p,'COORDS',{'C1','C2','C3'});
addParameter(p,'RANDOM_SLOPES',true);
parse(p,varargin{:});
cfg = p.Results;
C = cfg.COORDS; k = numel(C);

% -------------------- Build working table --------------------
S = load(cfg.MAT); T = S.tbl2;
W = table();
W.relief = -T.pain;                         % larger = more relief
W.StudyID = categorical(T.StudyID);   % numeric w/ gaps -> must be a grouping variable
Z = T{:,C};
Z = (Z - mean(Z)) ./ std(Z);                % global z-score -> standardized betas
for j = 1:k, W.(C{j}) = Z(:,j); end

fprintf('=== READOUT MAP w  (spec Step A, within-study LMM) ===\n');
fprintf('COORDS = {%s} | relief = -pain (standardized coords)\n\n', strjoin(C,','));

% -------------------- Model 1: random intercept --------------------
rhs = strjoin(C,' + ');
f1  = sprintf('relief ~ %s + (1|StudyID)', rhs);
lme = fitlme(W, f1);
fe  = lme.Coefficients;

fprintf('--- within-study readout (random intercept):  %s\n', f1);
w = zeros(k,1); se = zeros(k,1);
for j = 1:k
    r = strcmp(fe.Name, C{j});
    w(j) = fe.Estimate(r); se(j) = fe.SE(r);
    fprintf('  %-8s beta=%+.4f  SE=%.4f  t=%+.2f  p=%.2e\n', ...
        C{j}, w(j), se(j), fe.tStat(r), fe.pValue(r));
end
what = w / norm(w);
fprintf('  ||w|| = %.4f    w_hat = [%s]\n', norm(w), strjoin(compose('%+.3f',what'),', '));
fprintf('  R^2 (ordinary) = %.3f    R^2 (adjusted) = %.3f\n', ...
    lme.Rsquared.Ordinary, lme.Rsquared.Adjusted);

% low-rank read: share of ||w|| carried by the single dominant coordinate
[~,dom] = max(abs(w));
fprintf('  dominant axis = %s  (|beta| = %.1f%% of ||w||_1)\n', ...
    C{dom}, 100*abs(w(dom))/sum(abs(w)));

% -------------------- Null space of w --------------------
if k >= 2
    Nsp = null(w');                          % k x (k-1) orthonormal basis
    fprintf('\n--- null space of readout (iso-analgesia directions), dim %d ---\n', size(Nsp,2));
    for c = 1:size(Nsp,2)
        fprintf('  null_%d = [%s]  (w . null = %.1e)\n', c, ...
            strjoin(compose('%+.3f',Nsp(:,c)'),', '), w'*Nsp(:,c));
    end
    out.null = Nsp;
end

% -------------------- Model 2: random slopes (universality) --------------------
if cfg.RANDOM_SLOPES && k <= 3
    f2 = sprintf('relief ~ %s + (%s|StudyID)', rhs, rhs);
    fprintf('\n--- readout universality (random slopes): %s\n', f2);
    try
        lme_reml = fitlme(W, f1, 'FitMethod','REML');   % same FitMethod for LR test
        lme2 = fitlme(W, f2, 'FitMethod','REML');
        cmp = compare(lme_reml, lme2, 'CheckNesting', true);
        fprintf('  random-slope model LR test vs intercept-only: p = %.2e\n', cmp.pValue(2));
        fprintf('  (significant => readout w varies by study; a reportable finding)\n');
        % spread of by-study slopes
        [B,Bnames] = randomEffects(lme2);
        for j = 1:k
            sl = B(strcmp(Bnames.Name, C{j}));
            fprintf('    %-8s by-study slope SD = %.4f (fixed %+.4f)\n', C{j}, std(sl), w(j));
        end
        out.lme_rs = lme2;
    catch ME
        fprintf('  random-slope fit failed (%s); reporting intercept-only w.\n', ME.message);
    end
end

% -------------------- Within vs ecological (cross-study) slope --------------------
G = groupsummary(W, 'StudyID', 'mean', [{'relief'} C]);
Xe = G{:, compose('mean_%s', string(C))};
ye = G.mean_relief;
be = [ones(size(Xe,1),1) Xe] \ ye;           % OLS on study means
fprintf('\n--- ecological (study-mean) slope, for contrast with within-study w ---\n');
for j = 1:k
    fprintf('  %-8s within=%+.4f   ecological=%+.4f\n', C{j}, w(j), be(j+1));
end
fprintf('  (divergence => Simpson-type effect; within-study w is the correct readout)\n');

out.w = w; out.what = what; out.se = se; out.lme = lme; out.COORDS = C;
end
