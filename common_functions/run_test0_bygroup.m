function R = run_test0_bygroup(varargin)
% RUN_TEST0_BYGROUP  Test 0 broken out per intervention group.
%
% Same relief-banded decoding as run_test0/run_test0_loso, but instead of
% collapsing to one balanced accuracy per band, it tracks PER-CLASS recall
% and aggregates across bands. Gives, for each intervention group:
%
%   Acc     - mean per-class recall across relief bands (relief-matched)
%   Floor   - permutation null (labels shuffled within band)
%   Ceiling - pooled decoder over the full relief spectrum at matched N
%   p       - permutation p-value for that group
%
% MODE:
%   'cv'   - ordinary stratified k-fold CV (the optimistic ceiling pass)
%   'loso' - leave-one-study-out (the honest pass)
%
% Usage:
%   Rcv   = run_test0_bygroup('MODE','cv');
%   Rloso = run_test0_bygroup('MODE','loso');

p = inputParser;
addParameter(p,'MAT','datamat.mat');
addParameter(p,'MODE','cv');
addParameter(p,'COORDS',{'C1','C2','C3'});
addParameter(p,'N_BINS',7);
addParameter(p,'KFOLD',5);
addParameter(p,'N_PERM',[]);
addParameter(p,'N_POOL',[]);
addParameter(p,'MIN_CLASS_PER_BIN',5);
% Class frequencies here reflect which studies were collected, not a
% meaningful base rate. With empirical priors LDA just predicts Placebo
% (N=611) and per-group recall becomes a readout of class size, so default
% to uniform priors.
addParameter(p,'PRIOR','uniform');
addParameter(p,'RNG_SEED',7);
parse(p,varargin{:});
cfg = p.Results;
if isempty(cfg.N_PERM), cfg.N_PERM = ternary(strcmpi(cfg.MODE,'loso'),200,500); end
if isempty(cfg.N_POOL), cfg.N_POOL = ternary(strcmpi(cfg.MODE,'loso'),100,200); end
rng(cfg.RNG_SEED);

GRP = {'Cognitive','Conditioning','Control','Mindfulness', ...
       'Placebo','Placebo_C','Remifentanil','Social'};
nG = numel(GRP);

S = load(cfg.MAT); T = S.tbl2;
relief = -T.pain;
D = T{:,GRP}; [~,gi] = max(D,[],2);
g   = categorical(gi,1:nG,GRP);
sid = categorical(T.StudyID);          % numeric w/ gaps -> force ID semantics
X   = T{:,cfg.COORDS};
isLOSO = strcmpi(cfg.MODE,'loso');
if ~isLOSO
    X = (X - mean(X))./std(X);         % global z-score (train-only inside LOSO)
end

edges = quantile(relief, linspace(0,1,cfg.N_BINS+1));
edges(1) = -inf; edges(end) = inf;
bin = discretize(relief, edges);

fprintf('=== TEST 0 BY GROUP | MODE=%s | COORDS={%s} | %d bands ===\n', ...
    upper(cfg.MODE), strjoin(cfg.COORDS,','), cfg.N_BINS);

% ---------- observed ----------
obs = band_pass(X, g, sid, bin, GRP, cfg, isLOSO, false);

% ---------- permutation floor ----------
nullM = nan(cfg.N_PERM, nG);
for r = 1:cfg.N_PERM
    nullM(r,:) = band_pass(X, g, sid, bin, GRP, cfg, isLOSO, true);
end
floorV = mean(nullM,1,'omitnan');

% ---------- pooled ceiling at matched N (full relief spectrum) ----------
poolM = nan(cfg.N_POOL, nG);
for r = 1:cfg.N_POOL
    rec = nan(cfg.N_BINS, nG);
    for b = 1:cfg.N_BINS
        nb = sum(bin==b);
        idx = randsample(numel(relief), nb);      % drawn across ALL bands
        rec(b,:) = one_cell(X(idx,:), g(idx), sid(idx), GRP, cfg, isLOSO, false);
    end
    poolM(r,:) = mean(rec,1,'omitnan');
end
ceilV = mean(poolM,1,'omitnan');

% ---------- p-values ----------
pv = nan(nG,1);
for i = 1:nG
    if isnan(obs(i)), continue; end
    pv(i) = (1 + sum(nullM(:,i) >= obs(i))) / (cfg.N_PERM + 1);
end

nPer = arrayfun(@(i) sum(g==GRP{i}), 1:nG)';
R = table(GRP(:), obs(:), floorV(:), ceilV(:), pv, nPer, ...
    'VariableNames', {'Group','Acc','Floor','Ceiling','p','N'});
R.Properties.UserData = cfg;

fprintf('\n%-14s %7s %7s %7s %8s %6s\n','group','acc','floor','ceil','p','N');
for i = 1:nG
    fprintf('%-14s %7.3f %7.3f %7.3f %8.4f %6d\n', ...
        GRP{i}, R.Acc(i), R.Floor(i), R.Ceiling(i), R.p(i), R.N(i));
end
fprintf('\nmean acc=%.3f  mean floor=%.3f  mean ceiling=%.3f\n', ...
    mean(R.Acc,'omitnan'), mean(R.Floor,'omitnan'), mean(R.Ceiling,'omitnan'));
end

% ================= helpers =================
function v = band_pass(X, g, sid, bin, GRP, cfg, isLOSO, permute_labels)
% one full pass over all relief bands -> per-group recall averaged over bands
    nb = cfg.N_BINS; rec = nan(nb, numel(GRP));
    for b = 1:nb
        m = (bin==b);
        rec(b,:) = one_cell(X(m,:), g(m), sid(m), GRP, cfg, isLOSO, permute_labels);
    end
    v = mean(rec,1,'omitnan');
end

function rec = one_cell(Xb, gb, sb, GRP, cfg, isLOSO, permute_labels)
% per-class recall within a single cell (band, or matched-N subsample)
    rec = nan(1, numel(GRP));
    gb = removecats(categorical(gb));
    % drop classes too thin to be a target
    cls = categories(gb);
    keep = true(size(gb));
    for i = 1:numel(cls)
        if sum(gb==cls{i}) < cfg.MIN_CLASS_PER_BIN, keep(gb==cls{i}) = false; end
    end
    gb = removecats(gb(keep)); Xb = Xb(keep,:); sb = removecats(categorical(sb(keep)));
    if numel(categories(gb)) < 2 || numel(gb) < 2*cfg.KFOLD, return; end
    if permute_labels, gb = gb(randperm(numel(gb))); end

    if isLOSO
        [yhat, scored] = loso_oof(Xb, gb, sb, cfg.PRIOR);
    else
        yhat = cv_oof(Xb, gb, cfg.KFOLD, cfg.PRIOR); scored = true(size(gb));
    end
    ytrue = string(gb); yhat = string(yhat);
    for i = 1:numel(GRP)
        mi = scored & (ytrue==GRP{i});
        if any(mi), rec(i) = mean(yhat(mi)==GRP{i}); end
    end
end

function yhat = cv_oof(X, y, k, prior)
    y = removecats(categorical(y));
    c = cvpartition(y,'KFold',k); yhat = y;
    for f = 1:k
        tr = training(c,f); te = test(c,f);
        mdl = fitcdiscr(X(tr,:), y(tr), 'DiscrimType','pseudoLinear', 'Prior', prior);
        yhat(te) = predict(mdl, X(te,:));
    end
end

function [yhat, scored] = loso_oof(X, y, sid, prior)
% leave-one-study-out within a cell; train-only standardization per fold
    y = removecats(categorical(y)); sid = removecats(categorical(sid));
    studies = categories(sid);
    yhat = y; scored = false(numel(y),1); ytrue = string(y);
    for s = 1:numel(studies)
        te = (sid==studies{s}); tr = ~te;
        gtr = removecats(y(tr));
        if numel(categories(gtr)) < 2, continue; end
        mu = mean(X(tr,:),1); sg = std(X(tr,:),0,1); sg(sg==0)=1;
        mdl = fitcdiscr((X(tr,:)-mu)./sg, gtr, 'DiscrimType','pseudoLinear', 'Prior', prior);
        pred = predict(mdl, (X(te,:)-mu)./sg);
        teIdx = find(te); trCls = string(categories(gtr));
        for j = 1:numel(teIdx)
            k = teIdx(j);
            if ismember(ytrue(k), trCls)     % class must be learnable from training
                yhat(k) = pred(j); scored(k) = true;
            end
        end
    end
end

function out = ternary(c,a,b), if c, out=a; else, out=b; end, end
