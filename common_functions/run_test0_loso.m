function results = run_test0_loso(varargin)
% RUN_TEST0_LOSO  Test 0, Stage 2: the honest leave-one-study-out pass.
%
% Same within-band decoding as run_test0, but validated leave-one-study-out:
% a held-out study's subjects are classified using LDA trained ONLY on OTHER
% studies (train-only standardization), within the same relief band. Because
% study is perfectly confounded with intervention, this is the test that
% separates a generalizable intervention signature from a per-study batch
% signature:
%   survives LOSO  -> relief-matched intervention identity is real (mechanism)
%   collapses to chance -> the Stage-1 ceiling was study/scanner batch effects
%
% Floor = subject-level label permutation within band, re-run through the
% full LOSO loop (a generalizable-signal null).
%
% Usage:
%   results = run_test0_loso();
%   results = run_test0_loso('COORDS',{'NPS','SIIPS'});

% -------------------- Configuration --------------------
p = inputParser;
addParameter(p,'MAT','datamat.mat');
addParameter(p,'COORDS',{'C1','C2','C3'});
addParameter(p,'N_BINS',7);
addParameter(p,'N_PERM',200);
addParameter(p,'MIN_CLASS_PER_BIN',5);   % class must have >=this in band to be a target
addParameter(p,'RNG_SEED',7);
parse(p,varargin{:});
cfg = p.Results;
rng(cfg.RNG_SEED);

GRP = {'Cognitive','Conditioning','Control','Mindfulness', ...
       'Placebo','Placebo_C','Remifentanil','Social'};

% -------------------- Load & prepare --------------------
S = load(cfg.MAT); T = S.tbl2;
relief = -T.pain;
D = T{:,GRP}; [~,gi] = max(D,[],2);
g   = categorical(gi, 1:numel(GRP), GRP);
sid = categorical(T.StudyID);              % StudyID is numeric w/ gaps -> force ID semantics
Xraw = T{:,cfg.COORDS};                    % standardized train-only inside CV

edges = quantile(relief, linspace(0,1,cfg.N_BINS+1));
edges(1) = -inf; edges(end) = inf;
bin = discretize(relief, edges);

fprintf('=== TEST 0 (Stage 2: LEAVE-ONE-STUDY-OUT) ===\n');
fprintf('COORDS = {%s} | %d equal-count relief bands | seed %d | no harmonization\n', ...
    strjoin(cfg.COORDS,','), cfg.N_BINS, cfg.RNG_SEED);
fprintf('N = %d, %d groups, %d studies\n\n', height(T), numel(GRP), numel(categories(removecats(sid))));

nb = cfg.N_BINS;
BA_loso = nan(nb,1); BA_null_mu = nan(nb,1); BA_null_p95 = nan(nb,1);
pval = nan(nb,1); Nused = nan(nb,1); Nexcl = nan(nb,1);
nCls = nan(nb,1); reliefMid = nan(nb,1);

for b = 1:nb
    m = (bin==b);
    idx = find(m);
    reliefMid(b) = median(relief(idx));

    % observed LOSO balanced accuracy
    [BA_loso(b), Nused(b), Nexcl(b), nCls(b)] = ...
        loso_band(Xraw(idx,:), g(idx), sid(idx), cfg.MIN_CLASS_PER_BIN);

    % permutation floor: shuffle labels among band subjects, re-run LOSO
    null = nan(cfg.N_PERM,1);
    gi_b = g(idx);
    for r = 1:cfg.N_PERM
        gp = gi_b(randperm(numel(gi_b)));
        null(r) = loso_band(Xraw(idx,:), gp, sid(idx), cfg.MIN_CLASS_PER_BIN);
    end
    BA_null_mu(b) = mean(null,'omitnan');
    BA_null_p95(b) = prctile(null,95);
    pval(b) = (1 + sum(null >= BA_loso(b))) / (cfg.N_PERM + 1);

    fprintf(['band %2d | relief~%+.2f | used=%3d (excl %2d) C=%d | LOSO=%.3f  ' ...
             'null=%.3f (p95=%.3f, p=%.3f)  chance=%.3f\n'], b, reliefMid(b), ...
             Nused(b), Nexcl(b), nCls(b), BA_loso(b), BA_null_mu(b), ...
             BA_null_p95(b), pval(b), 1/max(nCls(b),1));
end

ok = ~isnan(BA_loso);
fprintf('\n--- aggregate over %d bands ---\n', sum(ok));
fprintf('mean LOSO BA        : %.3f\n', mean(BA_loso(ok)));
fprintf('mean permutation BA : %.3f\n', mean(BA_null_mu(ok)));
fprintf('LOSO - null         : %.3f\n', mean(BA_loso(ok))-mean(BA_null_mu(ok)));
fprintf(['Read: LOSO >> null -> generalizable relief-matched intervention signature; ' ...
         'LOSO ~= null -> Stage-1 ceiling was study batch.\n']);

results = table((1:nb)', reliefMid, Nused, Nexcl, nCls, BA_loso, BA_null_mu, ...
    BA_null_p95, pval, 'VariableNames', {'band','reliefMid','Nused','Nexcl', ...
    'nClasses','BA_loso','BA_null','BA_null_p95','p'});
results.Properties.UserData = cfg;
end

% ==================== helpers ====================
function [ba, nUsed, nExcl, nCls] = loso_band(X, g, sid, minN)
% Leave-one-study-out within a single relief band.
%  - target classes = groups with >= minN members in the band
%  - a held-out study's subject is scored only if its true class is present
%    in the training studies (else excluded, counted in nExcl)
%  - standardization fit on training subjects only, per held-out study
    g = removecats(categorical(g));
    sid = removecats(categorical(sid));

    % eligible target classes in this band
    cls = categories(g);
    tgt = cls(cellfun(@(c) sum(g==c) >= minN, cls));
    if numel(tgt) < 2
        ba = NaN; nUsed = 0; nExcl = sum(ismember(string(g),string(cls))); nCls = numel(tgt);
        return;
    end

    studies = categories(sid);
    yhat = strings(numel(g),1); ytrue = string(g); scored = false(numel(g),1);

    for s = 1:numel(studies)
        te = (sid==studies{s});
        tr = ~te;
        gtr = removecats(g(tr));
        trCls = categories(gtr);
        if numel(trCls) < 2, continue; end             % need >=2 classes to train

        % train-only standardization
        mu = mean(X(tr,:),1); sg = std(X(tr,:),0,1); sg(sg==0)=1;
        Xtr = (X(tr,:)-mu)./sg; Xte = (X(te,:)-mu)./sg;

        % uniform priors: class size here reflects study collection, not base rate
        mdl = fitcdiscr(Xtr, gtr, 'DiscrimType','pseudoLinear', 'Prior','uniform');
        pred = string(predict(mdl, Xte));

        teIdx = find(te);
        for j = 1:numel(teIdx)
            k = teIdx(j);
            % score only if this subject's true class is an eligible target
            % AND was available to be learned in training
            if ismember(ytrue(k), tgt) && ismember(ytrue(k), string(trCls))
                yhat(k) = pred(j); scored(k) = true;
            end
        end
    end

    % balanced accuracy over eligible target classes that were actually scored
    present = intersect(tgt, unique(ytrue(scored)), 'stable');
    rec = zeros(numel(present),1);
    for i = 1:numel(present)
        mi = scored & (ytrue==present{i});
        rec(i) = mean(yhat(mi)==present{i});
    end
    ba = mean(rec);
    nUsed = sum(scored); nExcl = sum(~scored); nCls = numel(present);
end
