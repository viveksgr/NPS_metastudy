function R = run_test1(varargin)
% RUN_TEST1  Spec §5 - variance-fraction / angle split of between-intervention
% separation into on-readout-axis vs null-space components.
%
% For every intervention pair (g,g') in whitened space:
%   Delta = zbar_g - zbar_g'
%   a = (w_hat' * Delta)^2                 on-axis  (1 dimension)
%   b = || Delta - (w_hat'*Delta) w_hat||^2 null-space (k-1 dimensions)
%   f = b/(a+b)                            null-space fraction
%   theta = angle between Delta and w_hat
%
% THREE BENCHMARKS (none of which are 0 / 50% / 90deg):
%   1. SOFT CLAIM - is there real null-space separation?
%        b vs its PERMUTATION null (shuffle group labels). b>0 is trivially
%        true for continuous data, so the null is the only meaningful floor.
%   2. HARD CLAIM - is the null space enriched relative to the readout axis?
%        Compare PER DIMENSION: r = [b/(k-1)] / a.  Isotropy => r = 1.
%        NB f = b/(a+b) > 0.5 is NOT evidence: isotropy already gives
%        f = (k-1)/k  (0.67 for k=3, 0.80 for k=5).
%   3. ANGLES - vs the isotropic angle distribution, simulated for this k.
%        For k=3 a random Delta sits at a median of 60deg from w_hat, not 0.
%
% MODE:
%   'all'         - all pairs over all studies. NB group means are also study
%                   means for single-arm studies, so Delta is batch-contaminated.
%   'withinstudy' - the 3 multi-arm studies only (5/7/9). Delta is computed
%                   within a study, so it is batch-free. 3 pairs, but clean.
%
% Usage:
%   R = run_test1();
%   R = run_test1('MODE','withinstudy');
%   R = run_test1('COORDS',{'NPS','SIIPS','C1','C2','C3'});

p = inputParser;
addParameter(p,'MAT','datamat.mat');
addParameter(p,'COORDS',{'C1','C2','C3'});
addParameter(p,'MODE','all');
addParameter(p,'N_PERM',2000);
addParameter(p,'N_BOOT',2000);
addParameter(p,'N_ISO',20000);
addParameter(p,'MIN_PER_GROUP',10);
addParameter(p,'RNG_SEED',7);
parse(p,varargin{:});
cfg = p.Results; rng(cfg.RNG_SEED);

GRP = {'Cognitive','Conditioning','Control','Mindfulness', ...
       'Placebo','Placebo_C','Remifentanil','Social'};
C = cfg.COORDS; k = numel(C);

S = load(cfg.MAT); T = S.tbl2;
relief = -T.pain;
D = T{:,GRP}; [~,gi] = max(D,[],2);
g   = categorical(gi,1:numel(GRP),GRP);
sid = categorical(T.StudyID);
X   = T{:,C}; X = (X-mean(X))./std(X);

fprintf('=== TEST 1: on-axis vs null-space separation ===\n');
fprintf('COORDS={%s} | k=%d (null dim=%d) | MODE=%s\n', strjoin(C,','), k, k-1, cfg.MODE);
fprintf('isotropic reference: f=(k-1)/k=%.3f, per-dim ratio r=1\n\n', (k-1)/k);

isWithin = strcmpi(cfg.MODE,'withinstudy');

% ---------------- geometry ----------------
if isWithin
    [P, pairName] = pairs_withinstudy(X, g, sid, relief, C, GRP, cfg);
else
    [P, pairName] = pairs_all(X, g, sid, relief, C, GRP, cfg, false);
end
if isempty(P), fprintf('no usable pairs\n'); R=[]; return; end

fprintf('%-34s %8s %8s %8s %8s\n','pair','a(on)','b(null)','f','angle');
for i = 1:numel(pairName)
    fprintf('%-34s %8.4f %8.4f %8.3f %7.1f\n', pairName{i}, ...
        P(i).a, P(i).b, P(i).f, P(i).ang);
end

fA = mean([P.a]); fB = mean([P.b]);
fMean = mean([P.f]); angMean = mean([P.ang]); angMed = median([P.ang]);
rPer  = (fB/(k-1)) / fA;

fprintf('\nmean a (on-axis)      = %.4f\n', fA);
fprintf('mean b (null space)   = %.4f\n', fB);
fprintf('mean f = b/(a+b)      = %.3f   [isotropy = %.3f]\n', fMean, (k-1)/k);
fprintf('per-dim ratio r       = %.2f    [isotropy = 1.00]\n', rPer);
fprintf('mean/median angle     = %.1f / %.1f deg\n', angMean, angMed);

% ---------------- benchmark 1: permutation null (soft claim) ----------------
permA = nan(cfg.N_PERM,1);
permB = nan(cfg.N_PERM,1); permF = nan(cfg.N_PERM,1); permR = nan(cfg.N_PERM,1);
% per-intervention nulls too, so each bar in the by-group figure can carry
% its own significance against the noise floor
permAg = nan(cfg.N_PERM, numel(GRP)); permBg = nan(cfg.N_PERM, numel(GRP));
for r = 1:cfg.N_PERM
    gp = g(randperm(numel(g)));
    if isWithin
        Pp = pairs_withinstudy(X, gp, sid, relief, C, GRP, cfg);
    else
        Pp = pairs_all(X, gp, sid, relief, C, GRP, cfg, true);
    end
    if isempty(Pp), continue; end
    permA(r) = mean([Pp.a]);
    permB(r) = mean([Pp.b]); permF(r) = mean([Pp.f]);
    permR(r) = (mean([Pp.b])/(k-1)) / mean([Pp.a]);
    [permAg(r,:), permBg(r,:)] = bygroup_means(Pp, GRP);
end
[obsAg, obsBg] = bygroup_means(P, GRP);
pB = (1+sum(permB >= fB))/(cfg.N_PERM+1);
pF = (1+sum(permF >= fMean))/(cfg.N_PERM+1);
pR = (1+sum(permR >= rPer))/(cfg.N_PERM+1);

% ---------------- benchmark 3: isotropic angle distribution ----------------
V = randn(cfg.N_ISO, k); V = V./vecnorm(V,2,2);
isoAng = acosd(abs(V(:,1)));                 % angle to an arbitrary fixed axis
pAng = mean(isoAng >= angMean);

fprintf('\n--- benchmark 1: permutation null (SOFT claim: real null-space separation?) ---\n');
fprintf('b: observed %.4f vs null %.4f (p=%.4f)  -> %s\n', fB, mean(permB,'omitnan'), pB, ...
    verdict(pB < 0.05, 'null-space separation EXCEEDS label noise', ...
                       'null-space separation NOT distinguishable from noise'));

fprintf('\n--- benchmark 2: isotropy (HARD claim: null space enriched per dimension?) ---\n');
fprintf('f: observed %.3f vs isotropic %.3f\n', fMean, (k-1)/k);
fprintf('r: observed %.2f vs isotropic 1.00 (perm p=%.4f) -> %s\n', rPer, pR, ...
    verdict(rPer > 1 && pR < 0.05, 'null space ENRICHED beyond dimension counting', ...
                                   'consistent with isotropy - no enrichment'));

% ---------------- noise floor: what permutation says a and b are ----------------
% Permuted Delta is pure sampling noise. After whitening, that noise is
% isotropic BY CONSTRUCTION - which is why the permutation null and the
% analytic isotropy reference land on the same place. The permuted values
% are also the noise contribution to the observed a and b, so subtracting
% them gives noise-corrected estimates.
mA = mean(permA,'omitnan'); mB = mean(permB,'omitnan');
fprintf('\n--- noise floor (permutation) vs analytic isotropy ---\n');
fprintf('perm mean a = %.4f | perm mean b = %.4f | perm b/a = %.2f  [isotropy predicts %.2f]\n', ...
    mA, mB, mB/mA, k-1);
fprintf('perm mean f = %.3f  [analytic isotropy %.3f]\n', mean(permF,'omitnan'), (k-1)/k);
fprintf('perm mean r = %.2f  [analytic isotropy 1.00]\n', mean(permR,'omitnan'));
fprintf('noise-corrected: a = %.4f, b = %.4f -> f = %.3f, r = %.2f\n', ...
    fA-mA, fB-mB, (fB-mB)/((fA-mA)+(fB-mB)), ((fB-mB)/(k-1))/(fA-mA));
fprintf('b as %% of a: raw %.0f%%  | per-dimension %.0f%%\n', 100*fB/fA, 100*rPer);

fprintf('\n--- benchmark 3: angles vs isotropic distribution for k=%d ---\n', k);
fprintf('observed mean angle %.1f deg | isotropic mean %.1f, median %.1f (p=%.4f)\n', ...
    angMean, mean(isoAng), median(isoAng), pAng);

% ---------------- bootstrap over studies ----------------
if ~isWithin
    bf = nan(cfg.N_BOOT,1); br = nan(cfg.N_BOOT,1);
    us = categories(removecats(sid));
    for bIdx = 1:cfg.N_BOOT
        pick = us(randi(numel(us), numel(us),1));
        m = false(numel(sid),1);
        for i = 1:numel(pick), m = m | (sid==pick{i}); end
        Pb = pairs_all(X(m,:), g(m), sid(m), relief(m), C, GRP, cfg, true);
        if isempty(Pb), continue; end
        bf(bIdx) = mean([Pb.f]);
        br(bIdx) = (mean([Pb.b])/(k-1))/mean([Pb.a]);
    end
    fprintf('\nbootstrap over studies: f = %.3f [%.3f, %.3f] | r = %.2f [%.2f, %.2f]\n', ...
        fMean, prctile(bf,2.5), prctile(bf,97.5), rPer, prctile(br,2.5), prctile(br,97.5));
end

% p-values against the permutation noise floor, overall and per intervention
pA_overall = (1+sum(permA >= fA))/(cfg.N_PERM+1);
pAg = nan(numel(GRP),1); pBg = nan(numel(GRP),1);
for i = 1:numel(GRP)
    if ~isnan(obsAg(i))
        pAg(i) = (1+sum(permAg(:,i) >= obsAg(i)))/(cfg.N_PERM+1);
        pBg(i) = (1+sum(permBg(:,i) >= obsBg(i)))/(cfg.N_PERM+1);
    end
end

R = struct('pairs',{P},'pairName',{pairName},'f',fMean,'r',rPer, ...
           'a',fA,'b',fB,'angle',angMean,'pB',pB,'pR',pR,'pF',pF, ...
           'pA',pA_overall,'obsAg',obsAg,'obsBg',obsBg,'pAg',pAg,'pBg',pBg, ...
           'permA',permA,'permB',permB,'GRP',{GRP},'cfg',cfg);
end

function [mA, mB] = bygroup_means(P, GRP)
% mean a and b over all pairs containing each intervention
    mA = nan(1,numel(GRP)); mB = nan(1,numel(GRP));
    if isempty(P), return; end
    g1 = {P.g1}; g2 = {P.g2}; aa = [P.a]; bb = [P.b];
    for i = 1:numel(GRP)
        m = strcmp(g1,GRP{i}) | strcmp(g2,GRP{i});
        if any(m), mA(i) = mean(aa(m)); mB(i) = mean(bb(m)); end
    end
end

% ==================== helpers ====================
function [P, nm] = pairs_all(X, g, sid, relief, C, GRP, cfg, quiet)
    P = struct('a',{},'b',{},'f',{},'ang',{},'g1',{},'g2',{}); nm = {};
    g = categorical(g);
    Sig = resid_cov(X, g, sid); Wm = inv_sqrtm(Sig); Z = X*Wm;
    what = fit_readout(Z, relief, sid, C, quiet);
    if isempty(what), return; end
    present = GRP(cellfun(@(c) sum(g==c) >= cfg.MIN_PER_GROUP, GRP));
    for i = 1:numel(present)
        for j = i+1:numel(present)
            d = (mean(Z(g==present{i},:),1) - mean(Z(g==present{j},:),1))';
            [P(end+1), s] = split_delta(d, what, present{i}, present{j}); %#ok<AGROW>
            nm{end+1} = s; %#ok<AGROW>
        end
    end
end

function [P, nm] = pairs_withinstudy(X, g, sid, relief, C, GRP, cfg)
% Delta computed WITHIN each multi-arm study -> batch-free.
% w_hat taken from all OTHER studies (w is universal, LR p=0.71), so the
% readout axis is not fitted on the same subjects the contrast uses.
    P = struct('a',{},'b',{},'f',{},'ang',{},'g1',{},'g2',{}); nm = {};
    g = categorical(g); sid = removecats(categorical(sid));
    studies = categories(sid);
    for s = 1:numel(studies)
        m = sid==studies{s};
        gg = removecats(g(m));
        cls = categories(gg);
        if numel(cls) < 2, continue; end
        if any(cellfun(@(c) sum(gg==c) < cfg.MIN_PER_GROUP, cls)), continue; end

        Xs = X(m,:); Xs = (Xs-mean(Xs))./std(Xs);
        Sig = resid_cov(Xs, gg, sid(m)); Wm = inv_sqrtm(Sig); Zs = Xs*Wm;

        % readout from the other studies
        o = ~m;
        Xo = X(o,:); Xo = (Xo-mean(Xo))./std(Xo);
        So = resid_cov(Xo, g(o), sid(o)); Zo = Xo*inv_sqrtm(So);
        what = fit_readout(Zo, relief(o), sid(o), C, true);
        if isempty(what), continue; end

        d = (mean(Zs(gg==cls{1},:),1) - mean(Zs(gg==cls{2},:),1))';
        [P(end+1), nmS] = split_delta(d, what, cls{1}, cls{2}); %#ok<AGROW>
        nm{end+1} = sprintf('study %s: %s', studies{s}, nmS); %#ok<AGROW>
    end
end

function [out, nm] = split_delta(d, what, g1, g2)
    on  = (what'*d);
    out.a = on^2;
    out.b = norm(d - on*what)^2;
    out.f = out.b/(out.a+out.b);
    out.ang = acosd(min(1,abs(on)/norm(d)));
    out.g1 = g1;                 % kept so plots can aggregate by intervention
    out.g2 = g2;
    nm = sprintf('%s vs %s', g1, g2);
end

function Sig = resid_cov(X, g, sid)
    g = removecats(categorical(g)); sid = removecats(categorical(sid));
    Rz = X; cellid = strcat(string(g),'_',string(sid)); u = unique(cellid);
    for i = 1:numel(u)
        m = (cellid==u(i));
        if sum(m) > 1, Rz(m,:) = X(m,:)-mean(X(m,:),1); else, Rz(m,:) = 0; end
    end
    Sig = cov(Rz) + 1e-8*eye(size(X,2));
end

function Wm = inv_sqrtm(Sig)
    [V,Dg] = eig((Sig+Sig')/2); d = max(diag(Dg),1e-10);
    Wm = V*diag(1./sqrt(d))*V';
end

function what = fit_readout(Z, y, sid, C, quiet)
    what = [];
    W = table(); W.relief = y; W.StudyID = removecats(categorical(sid));
    for j = 1:size(Z,2), W.(C{j}) = Z(:,j); end
    try
        if numel(categories(W.StudyID)) > 1
            lme = fitlme(W, sprintf('relief ~ %s + (1|StudyID)', strjoin(C,' + ')));
        else
            lme = fitlme(W, sprintf('relief ~ %s', strjoin(C,' + ')));
        end
        fe = lme.Coefficients; w = zeros(size(Z,2),1);
        for j = 1:size(Z,2), w(j) = fe.Estimate(strcmp(fe.Name,C{j})); end
        what = w/norm(w);
    catch ME
        if ~quiet, fprintf('  readout fit failed: %s\n', ME.message); end
    end
end

function s = verdict(c,a,b), if c, s=a; else, s=b; end, end
