function results = run_test0_withinstudy(varargin)
% RUN_TEST0_WITHINSTUDY  Test 0, de-confounded: within-study, relief-matched
% arm decoding, using only the studies that contain >1 intervention arm.
%
% Within one study both arms share scanner, paradigm, preprocessing and
% population, so any arm separation is batch-free by construction. No
% cross-study generalization and no harmonization are involved.
%
% Two corrections relative to the first version of this script:
%   1. PAIRED SUBJECTS. These are crossover designs - arms are index-aligned,
%      so the k-th row of arm A and the k-th row of arm B are the SAME person.
%      Ordinary k-fold CV leaks: one of a subject's rows can sit in training
%      while the other is tested, letting the classifier recognize the person
%      rather than the intervention. We therefore use LEAVE-ONE-SUBJECT-OUT
%      CV, which holds out both of a subject's rows together. This also
%      removes the small-N fold instability of 5-fold at N=14-15.
%   2. UNIFORM PRIORS (see cv helper).
%
% Usage:
%   results = run_test0_withinstudy();
%   results = run_test0_withinstudy('COORDS',{'NPS','SIIPS'});

p = inputParser;
addParameter(p,'MAT','datamat.mat');
addParameter(p,'COORDS',{'C1','C2','C3'});
addParameter(p,'N_BINS',3);           % few bands: per-study N is small
addParameter(p,'N_PERM',1000);
addParameter(p,'N_POOL',500);
addParameter(p,'PRIOR','uniform');
addParameter(p,'RNG_SEED',7);
parse(p,varargin{:});
cfg = p.Results; rng(cfg.RNG_SEED);

GRP = {'Cognitive','Conditioning','Control','Mindfulness', ...
       'Placebo','Placebo_C','Remifentanil','Social'};

S = load(cfg.MAT); T = S.tbl2;
relief = -T.pain;
D = T{:,GRP}; [~,gi] = max(D,[],2);
g = categorical(gi,1:numel(GRP),GRP);
sid = categorical(T.StudyID);

fprintf('=== TEST 0, WITHIN-STUDY (de-confounded, paired subjects) ===\n');
fprintf('COORDS={%s} | %d bands/study | leave-one-subject-out CV | %s priors | NO harmonization\n\n', ...
    strjoin(cfg.COORDS,','), cfg.N_BINS, cfg.PRIOR);

studies = categories(removecats(sid));
rows = {};
for i = 1:numel(studies)
    m = sid==studies{i};
    gg = removecats(g(m));
    if numel(categories(gg)) < 2, continue; end

    X = T{m,cfg.COORDS}; X = (X - mean(X))./std(X);   % standardize within study
    rel = relief(m); arms = gg;
    armNames = categories(arms);
    subj = pair_subjects(arms);                        % index-aligned crossover
    nSubj = numel(unique(subj));
    fprintf('--- study %s | N=%d rows, %d subjects | %s vs %s ---\n', ...
        studies{i}, sum(m), nSubj, armNames{1}, armNames{2});

    e = quantile(rel, linspace(0,1,cfg.N_BINS+1)); e(1)=-inf; e(end)=inf;
    bin = discretize(rel, e);

    for b = 1:cfg.N_BINS
        mb = (bin==b);
        Xb = X(mb,:); ab = removecats(arms(mb)); sb = subj(mb);
        if numel(categories(ab)) < 2 || numel(unique(sb)) < 4
            fprintf('   band %d: N=%d too small - skipped\n', b, sum(mb)); continue;
        end
        ba = loso_subj_balacc(Xb, ab, sb, cfg.PRIOR);

        null = nan(cfg.N_PERM,1);
        for r = 1:cfg.N_PERM
            null(r) = loso_subj_balacc(Xb, ab(randperm(numel(ab))), sb, cfg.PRIOR);
        end
        pv = (1+sum(null>=ba))/(cfg.N_PERM+1);

        % pooled within-study reference at matched N (full relief spectrum)
        pool = nan(cfg.N_POOL,1);
        nb = sum(mb);
        for r = 1:cfg.N_POOL
            idx = randsample(numel(rel), nb);
            ap = removecats(arms(idx));
            if numel(categories(ap))>=2 && numel(unique(subj(idx)))>=4
                pool(r) = loso_subj_balacc(X(idx,:), ap, subj(idx), cfg.PRIOR);
            end
        end
        bp = mean(pool,'omitnan');

        fprintf(['   band %d | relief~%+.2f | N=%3d rows/%2d subj | within=%.3f  ' ...
                 'null=%.3f (p=%.4f)  pooled=%.3f  chance=0.500\n'], ...
                 b, median(rel(mb)), nb, numel(unique(sb)), ba, mean(null), pv, bp);
        rows(end+1,:) = {studies{i}, sprintf('%s_vs_%s',armNames{1},armNames{2}), ...
                         b, median(rel(mb)), nb, ba, mean(null), pv, bp}; %#ok<AGROW>
    end

    ba_all = loso_subj_balacc(X, arms, subj, cfg.PRIOR);
    fprintf('   [whole study, unbanded] BA=%.3f\n\n', ba_all);
end

results = cell2table(rows, 'VariableNames', {'study','contrast','band','reliefMid', ...
    'N','BA_within','BA_null','p','BA_pooled'});

if ~isempty(rows)
    fprintf('=== aggregate over %d relief-matched within-study bands ===\n', height(results));
    fprintf('mean within-band BA : %.3f\n', mean(results.BA_within));
    fprintf('mean permutation BA : %.3f\n', mean(results.BA_null));
    fprintf('mean pooled BA      : %.3f\n', mean(results.BA_pooled));
    fprintf('within - null       : %+.3f\n', mean(results.BA_within)-mean(results.BA_null));
    fprintf('within - pooled     : %+.3f\n', mean(results.BA_within)-mean(results.BA_pooled));
end
end

% ================= helpers =================
function subj = pair_subjects(arms)
% Crossover designs with index-aligned arms: the k-th subject of arm A is the
% same person as the k-th subject of arm B. Number subjects by position
% within each arm so both of a person's rows share an id.
    subj = nan(numel(arms),1);
    cls = categories(removecats(arms));
    for i = 1:numel(cls)
        idx = find(arms==cls{i});
        subj(idx) = 1:numel(idx);
    end
end

function ba = loso_subj_balacc(X, y, subj, prior)
% Leave-one-SUBJECT-out CV: both rows of a held-out subject go to test, so a
% person is never split across train/test. Uniform priors so the arm with
% more rows in a band cannot be favoured.
    y = removecats(categorical(y)); cls = categories(y);
    us = unique(subj); yhat = y; scored = false(numel(y),1);
    for i = 1:numel(us)
        te = (subj==us(i)); tr = ~te;
        ytr = removecats(y(tr));
        if numel(categories(ytr)) < 2, continue; end
        mdl = fitcdiscr(X(tr,:), ytr, 'DiscrimType','pseudoLinear', 'Prior', prior);
        yhat(te) = predict(mdl, X(te,:));
        scored(te) = true;
    end
    rec = nan(numel(cls),1);
    for i = 1:numel(cls)
        mi = scored & (y==cls{i});
        if any(mi), rec(i) = mean(yhat(mi)==cls{i}); end
    end
    ba = mean(rec,'omitnan');
end
