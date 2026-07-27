function results = run_test0(varargin)
% RUN_TEST0  Test 0 (model-free degeneracy), Stage 1: the "ceiling" pass.
%
% Bins subjects into narrow bands of ACTUAL relief (y = -pain), then within
% each band decodes intervention group from neural coordinates. Compares:
%   (1) within-band balanced accuracy (relief-matched decoding)
%   (2) a permutation FLOOR (labels shuffled within band)
%   (3) a pooled CEILING trained across the full pain spectrum, subsampled
%       to the per-band N so the gap reflects lost relief info, not lost N.
%
% This stage uses ordinary stratified k-fold CV (the achievable ceiling).
% Stage 2 (LOSO) will re-run the honest, study-out version.
%
% Usage:
%   results = run_test0();                          % defaults, COORDS={C1,C2,C3}
%   results = run_test0('COORDS',{'NPS','SIIPS'});  % robustness set
%
% ---- Data contract: datamat.mat -> table tbl2 (1645 x 16) ----
%   StudyID (39), StudyLabel (8), 8 one-hot intervention dummies,
%   pain (=> relief = -pain), NPS SIIPS C1 C2 C3.
%   NOTE: study is perfectly confounded with intervention (each study = one
%   group). Within-band decoding here therefore also reflects study/scanner
%   identity; that is what Stage-2 LOSO is designed to expose.

% -------------------- Configuration --------------------
p = inputParser;
addParameter(p,'MAT','datamat.mat');
addParameter(p,'COORDS',{'C1','C2','C3'});      % primary set
addParameter(p,'N_BINS',7);                     % equal-count relief bands
addParameter(p,'KFOLD',5);
addParameter(p,'N_PERM',1000);                  % permutation-floor draws
addParameter(p,'N_POOL_REP',200);              % pooled-subsample repeats
addParameter(p,'MIN_CLASS_PER_BIN',5);          % drop classes thinner than this in a band
addParameter(p,'RNG_SEED',7);
parse(p,varargin{:});
cfg = p.Results;
rng(cfg.RNG_SEED);

GRP = {'Cognitive','Conditioning','Control','Mindfulness', ...
       'Placebo','Placebo_C','Remifentanil','Social'};

% -------------------- Load & prepare --------------------
S = load(cfg.MAT); T = S.tbl2;
relief = -T.pain;                               % flip: larger = more relief
D = T{:,GRP};                                    % one-hot
[~,gi] = max(D,[],2);
g = categorical(gi, 1:numel(GRP), GRP);          % intervention label per subject

X = T{:,cfg.COORDS};
X = (X - mean(X)) ./ std(X);                      % global z-score (ceiling pass)

% Equal-count (quantile) relief bins
edges = quantile(relief, linspace(0,1,cfg.N_BINS+1));
edges(1) = -inf; edges(end) = inf;               % capture extremes
bin = discretize(relief, edges);

fprintf('=== TEST 0 (ceiling pass) ===\n');
fprintf('COORDS = {%s} | %d equal-count relief bands | %d-fold CV | seed %d\n', ...
    strjoin(cfg.COORDS,','), cfg.N_BINS, cfg.KFOLD, cfg.RNG_SEED);
fprintf('N = %d subjects, %d groups, %d studies\n\n', height(T), numel(GRP), numel(unique(T.StudyID)));

% -------------------- Per-band decoding --------------------
nb = cfg.N_BINS;
BA_within = nan(nb,1); BA_pool = nan(nb,1);
BA_null_mu = nan(nb,1); BA_null_p95 = nan(nb,1); pval = nan(nb,1);
Nband = nan(nb,1); nCls = nan(nb,1); reliefMid = nan(nb,1);

for b = 1:nb
    m = (bin==b);
    Xb = X(m,:); gb = removecats(g(m));
    % keep only classes with enough members in this band
    [gb, Xb] = drop_thin_classes(gb, Xb, cfg.MIN_CLASS_PER_BIN);
    Nband(b) = numel(gb); nCls(b) = numel(categories(gb));
    reliefMid(b) = median(relief(m));
    if nCls(b) < 2 || Nband(b) < 2*cfg.KFOLD
        fprintf('band %d: too few usable subjects/classes (N=%d, C=%d) - skipped\n', b, Nband(b), nCls(b));
        continue;
    end

    % (1) within-band balanced accuracy
    BA_within(b) = cv_balacc(Xb, gb, cfg.KFOLD);

    % (2) permutation floor (shuffle labels within band)
    null = nan(cfg.N_PERM,1);
    for r = 1:cfg.N_PERM
        null(r) = cv_balacc(Xb, gb(randperm(numel(gb))), cfg.KFOLD);
    end
    BA_null_mu(b) = mean(null); BA_null_p95(b) = prctile(null,95);
    pval(b) = (1 + sum(null >= BA_within(b))) / (cfg.N_PERM + 1);

    % (3) pooled ceiling, subsampled to this band's N (full spectrum)
    pool = nan(cfg.N_POOL_REP,1);
    for r = 1:cfg.N_POOL_REP
        idx = randsample(height(T), Nband(b));
        gp = removecats(g(idx)); Xp = X(idx,:);
        [gp, Xp] = drop_thin_classes(gp, Xp, cfg.MIN_CLASS_PER_BIN);
        if numel(categories(gp)) >= 2 && numel(gp) >= 2*cfg.KFOLD
            pool(r) = cv_balacc(Xp, gp, cfg.KFOLD);
        end
    end
    BA_pool(b) = mean(pool,'omitnan');

    fprintf(['band %2d | relief~%+.2f | N=%3d C=%d | within=%.3f  null=%.3f (p95=%.3f, p=%.3f)  ' ...
             'pooled=%.3f  chance=%.3f\n'], b, reliefMid(b), Nband(b), nCls(b), ...
             BA_within(b), BA_null_mu(b), BA_null_p95(b), pval(b), BA_pool(b), 1/nCls(b));
end

% -------------------- Aggregate --------------------
ok = ~isnan(BA_within);
fprintf('\n--- aggregate over %d usable bands ---\n', sum(ok));
fprintf('mean within-band BA : %.3f\n', mean(BA_within(ok)));
fprintf('mean permutation BA : %.3f\n', mean(BA_null_mu(ok)));
fprintf('mean pooled BA      : %.3f\n', mean(BA_pool(ok)));
fprintf('within-vs-pooled gap: %.3f  (within/pooled = %.2f)\n', ...
    mean(BA_pool(ok))-mean(BA_within(ok)), mean(BA_within(ok))/mean(BA_pool(ok)));
fprintf(['Read: within~=pooled -> relief-invariant identity; ' ...
         'within>>chance but <<pooled -> some decodability rode on relief.\n']);

results = table((1:nb)', reliefMid, Nband, nCls, BA_within, BA_null_mu, ...
    BA_null_p95, pval, BA_pool, ...
    'VariableNames', {'band','reliefMid','N','nClasses','BA_within', ...
    'BA_null','BA_null_p95','p','BA_pooled'});
results.Properties.UserData = cfg;
end

% ==================== helpers ====================
function ba = cv_balacc(X, y, k)
% Stratified k-fold, pooled-covariance (pseudo)linear LDA, balanced accuracy
% over out-of-fold predictions. pseudoLinear tolerates small/singular bands.
% Uniform priors: class sizes here reflect which studies were collected, not
% a base rate. With empirical priors LDA just predicts Placebo (N=611).
    y = removecats(categorical(y));
    cls = categories(y);
    c = cvpartition(y, 'KFold', k);
    yhat = y;                                   % preallocate as categorical
    for f = 1:k
        tr = training(c,f); te = test(c,f);
        mdl = fitcdiscr(X(tr,:), y(tr), 'DiscrimType','pseudoLinear', 'Prior','uniform');
        yhat(te) = predict(mdl, X(te,:));
    end
    rec = zeros(numel(cls),1);
    for i = 1:numel(cls)
        mi = (y==cls{i});
        rec(i) = mean(yhat(mi)==cls{i});
    end
    ba = mean(rec);
end

function [y2, X2] = drop_thin_classes(y, X, minN)
% Remove classes with fewer than minN members from the decoding set.
    y = removecats(categorical(y));
    cls = categories(y);
    keep = true(size(y));
    for i = 1:numel(cls)
        if sum(y==cls{i}) < minN
            keep(y==cls{i}) = false;
        end
    end
    y2 = removecats(y(keep)); X2 = X(keep,:);
end
