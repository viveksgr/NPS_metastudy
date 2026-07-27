function R = run_test2(varargin)
% RUN_TEST2  Spec §6 - decode intervention from the analgesia null space.
%
% Per subject, in whitened space:
%   a_i      = w_hat' * z_i            (readout coordinate)
%   z_i^perp = z_i - a_i * w_hat       (null-space residual)
% and decode group from three feature sets with the SAME estimator and folds:
%   A_full    - full z                    (3 dims)
%   A_readout - a alone                   (1 dim)
%   A_null    - null-space coords alone   (2 dims)
%
% Prediction (spec): A_null ~= A_full >> A_readout.
%
% DIMENSIONALITY CONTROL. A_null (2d) beating A_readout (1d) is partly
% arithmetic. We therefore also decode from random 1-D and random 2-D
% projections of the same whitened space:
%   A_readout ~= A_rand1 -> readout axis is not special, 1-D is just weak
%   A_readout <<  A_rand1 -> readout axis is SPECIFICALLY uninformative
% Only the second supports the spec's claim.
%
% IMPLEMENTATION NOTES
%   - Sigma_resid, standardization and w are fit on TRAINING studies only
%     (spec §8): noise in w leaks on-axis signal into the null projection.
%   - The null space is represented by an explicit orthonormal basis, so
%     [a, null_coords] is an orthogonal rotation of z and full = readout +
%     null exactly. (LDA is affine-invariant, so A_full is the same whether
%     computed on z or on the rotated coordinates.)
%   - Uniform priors + balanced accuracy + permutation floor, as in Test 0.
%
% Usage:
%   R = run_test2();                            % LOSO over all studies, 8-way
%   R = run_test2('COORDS',{'NPS','SIIPS'});    % robustness set (null dim = 1)

p = inputParser;
addParameter(p,'MAT','datamat.mat');
addParameter(p,'COORDS',{'C1','C2','C3'});
addParameter(p,'N_PERM',500);
addParameter(p,'N_RAND',200);
addParameter(p,'N_BOOT',1000);
addParameter(p,'PRIOR','uniform');
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
X   = T{:,C};

fprintf('=== TEST 2: decode from the null space ===\n');
fprintf('COORDS={%s} (k=%d, null dim=%d) | LOSO | %s priors | seed %d\n\n', ...
    strjoin(C,','), k, k-1, cfg.PRIOR, cfg.RNG_SEED);

studies = categories(removecats(sid));
nS = numel(studies);

% ---------------- pass 1: per-fold geometry (cached) ----------------
F = struct('tr',{},'te',{},'Ztr',{},'Zte',{},'what',{},'Nb',{});
wAll = nan(nS,k);
for s = 1:nS
    te = (sid==studies{s}); tr = ~te;
    if numel(categories(removecats(g(tr)))) < 2, continue; end

    % standardize on training studies only
    mu = mean(X(tr,:),1); sg = std(X(tr,:),0,1); sg(sg==0)=1;
    Xtr = (X(tr,:)-mu)./sg;  Xte = (X(te,:)-mu)./sg;

    % Sigma_resid: covariance after removing group x study cell means (train only)
    Sig = resid_cov(Xtr, g(tr), sid(tr));
    Wm  = inv_sqrtm(Sig);
    Ztr = Xtr*Wm;  Zte = Xte*Wm;

    % readout w in whitened space, train only, study random intercept
    what = fit_readout(Ztr, relief(tr), sid(tr), C);
    Nb   = null(what');                       % k x (k-1) orthonormal

    F(end+1) = struct('tr',tr,'te',te,'Ztr',Ztr,'Zte',Zte,'what',what,'Nb',Nb); %#ok<AGROW>
    wAll(s,:) = what';
end
fprintf('fitted %d LOSO folds. mean w_hat = [%s]\n', numel(F), ...
    strjoin(compose('%+.3f', mean(wAll,1,'omitnan')),', '));
fprintf('w_hat stability across folds (SD): [%s]\n\n', ...
    strjoin(compose('%.3f', std(wAll,0,1,'omitnan')),', '));

% ---------------- pass 2: decode from each feature set ----------------
sets = {'full','readout','null'};
acc = struct(); oof = struct();
accG = nan(numel(GRP), numel(sets));            % per-intervention recall
for i = 1:numel(sets)
    [acc.(sets{i}), oof.(sets{i}), rg] = loso_decode(F, g, sets{i}, [], cfg.PRIOR, GRP);
    accG(:,i) = rg(:);
end

% random projections (dimensionality controls)
r1 = nan(cfg.N_RAND,1); r2 = nan(cfg.N_RAND,1);
for r = 1:cfg.N_RAND
    r1(r) = loso_decode(F, g, 'rand', 1,   cfg.PRIOR, GRP);
    r2(r) = loso_decode(F, g, 'rand', k-1, cfg.PRIOR, GRP);
end

% permutation floors
perm = struct(); permG = struct();
for i = 1:numel(sets)
    v = nan(cfg.N_PERM,1); vG = nan(cfg.N_PERM, numel(GRP));
    for r = 1:cfg.N_PERM
        [v(r), ~, rg] = loso_decode(F, g(randperm(numel(g))), sets{i}, [], cfg.PRIOR, GRP);
        vG(r,:) = rg(:)';
    end
    perm.(sets{i}) = v; permG.(sets{i}) = vG;
end

% p-values against the permutation floor, per intervention and per feature set
pG = nan(numel(GRP), numel(sets));
for i = 1:numel(sets)
    for j = 1:numel(GRP)
        if ~isnan(accG(j,i))
            pG(j,i) = (1+sum(permG.(sets{i})(:,j) >= accG(j,i)))/(cfg.N_PERM+1);
        end
    end
end

% ---------------- bootstrap CIs over studies ----------------
bs = struct(); bsG = struct();
for i = 1:numel(sets)
    [bs.(sets{i}), bsG.(sets{i})] = boot_over_studies(oof.(sets{i}), sid, GRP, cfg.N_BOOT);
end

% ---------------- report ----------------
fprintf('%-12s %7s %7s %16s %10s\n','features','dims','bal.acc','95%% CI','perm p');
for i = 1:numel(sets)
    nm = sets{i};
    d  = pick(nm, k);
    pv = (1+sum(perm.(nm) >= acc.(nm)))/(cfg.N_PERM+1);
    ci = prctile(bs.(nm),[2.5 97.5]);
    fprintf('%-12s %7d %7.3f   [%.3f, %.3f] %10.4f\n', nm, d, acc.(nm), ci(1), ci(2), pv);
end
fprintf('%-12s %7d %7.3f   [%.3f, %.3f]\n','rand-1d',1,mean(r1),prctile(r1,2.5),prctile(r1,97.5));
fprintf('%-12s %7d %7.3f   [%.3f, %.3f]\n',sprintf('rand-%dd',k-1),k-1,mean(r2), ...
    prctile(r2,2.5),prctile(r2,97.5));
fprintf('%-12s %7s %7.3f\n','perm floor','-',mean(perm.full));

fprintf('\n--- readings ---\n');
fprintf('A_null / A_full     = %.2f   (spec predicts ~1)\n', acc.null/acc.full);
fprintf('A_null / A_readout  = %.2f   (spec predicts >>1)\n', acc.null/acc.readout);
fprintf('A_readout vs rand-1d: %.3f vs %.3f  -> %s\n', acc.readout, mean(r1), ...
    ternary(acc.readout < prctile(r1,5), ...
    'readout axis SPECIFICALLY uninformative (below random 1-D)', ...
    'readout axis NOT special; 1-D is simply weak'));
fprintf('A_null vs rand-%dd   : %.3f vs %.3f\n', k-1, acc.null, mean(r2));

R = struct('acc',acc,'perm',perm,'rand1',r1,'rand2',r2,'boot',bs, ...
           'accG',accG,'pG',pG,'permG',permG,'bootG',bsG,'sets',{sets}, ...
           'GRP',{GRP},'w',wAll,'cfg',cfg);
end

% ==================== helpers ====================
function Sig = resid_cov(X, g, sid)
% covariance of X after removing group x study cell means (spec §2.3)
    g = removecats(categorical(g)); sid = removecats(categorical(sid));
    Rz = X;
    cellid = strcat(string(g),'_',string(sid));
    u = unique(cellid);
    for i = 1:numel(u)
        m = (cellid==u(i));
        if sum(m) > 1, Rz(m,:) = X(m,:) - mean(X(m,:),1); else, Rz(m,:) = 0; end
    end
    Sig = cov(Rz);
    Sig = Sig + 1e-8*eye(size(Sig));      % guard
end

function Wm = inv_sqrtm(Sig)
    [V,Dg] = eig((Sig+Sig')/2);
    d = max(diag(Dg), 1e-10);
    Wm = V*diag(1./sqrt(d))*V';
end

function what = fit_readout(Z, y, sid, C)
% within-study readout in whitened space, study random intercept
    W = table(); W.relief = y; W.StudyID = removecats(categorical(sid));
    for j = 1:size(Z,2), W.(C{j}) = Z(:,j); end
    lme = fitlme(W, sprintf('relief ~ %s + (1|StudyID)', strjoin(C,' + ')));
    fe = lme.Coefficients;
    w = zeros(size(Z,2),1);
    for j = 1:size(Z,2), w(j) = fe.Estimate(strcmp(fe.Name,C{j})); end
    what = w/norm(w);
end

function [ba, oofOut, rec] = loso_decode(F, g, mode, dim, prior, GRP)
% decode group across cached LOSO folds using the requested feature set
    n = numel(g); yhat = strings(n,1); scored = false(n,1);
    for f = 1:numel(F)
        [Ftr, Fte] = features(F(f), mode, dim);
        gtr = removecats(g(F(f).tr));
        if numel(categories(gtr)) < 2, continue; end
        mdl = fitcdiscr(Ftr, gtr, 'DiscrimType','pseudoLinear','Prior',prior);
        pred = string(predict(mdl, Fte));
        idx = find(F(f).te); trCls = string(categories(gtr));
        ok = ismember(string(g(idx)), trCls);
        yhat(idx(ok)) = pred(ok); scored(idx(ok)) = true;
    end
    ytrue = string(g);
    rec = nan(numel(GRP),1);
    for i = 1:numel(GRP)
        m = scored & (ytrue==GRP{i});
        if any(m), rec(i) = mean(yhat(m)==GRP{i}); end
    end
    ba = mean(rec,'omitnan');
    if nargout > 1, oofOut = struct('yhat',yhat,'ytrue',ytrue,'scored',scored); end
end

function [Ftr, Fte] = features(Fs, mode, dim)
    switch mode
        case 'full'                       % orthogonal rotation of z
            Ftr = [Fs.Ztr*Fs.what, Fs.Ztr*Fs.Nb];
            Fte = [Fs.Zte*Fs.what, Fs.Zte*Fs.Nb];
        case 'readout'
            Ftr = Fs.Ztr*Fs.what;  Fte = Fs.Zte*Fs.what;
        case 'null'
            Ftr = Fs.Ztr*Fs.Nb;    Fte = Fs.Zte*Fs.Nb;
        case 'rand'
            P = orth(randn(size(Fs.Ztr,2), dim));
            Ftr = Fs.Ztr*P;        Fte = Fs.Zte*P;
    end
end

function [v, vG] = boot_over_studies(oof, sid, GRP, nboot)
% resample STUDIES with replacement; recompute balanced accuracy from the
% cached out-of-fold predictions (study is the resampling unit, spec §8)
    sid = removecats(categorical(sid)); us = categories(sid);
    v = nan(nboot,1); vG = nan(nboot, numel(GRP));
    for b = 1:nboot
        pick_ = us(randi(numel(us), numel(us),1));
        m = false(numel(sid),1);
        for i = 1:numel(pick_), m = m | (sid==pick_{i}); end
        m = m & oof.scored;
        rec = nan(numel(GRP),1);
        for i = 1:numel(GRP)
            mi = m & (oof.ytrue==GRP{i});
            if any(mi), rec(i) = mean(oof.yhat(mi)==GRP{i}); end
        end
        vG(b,:) = rec(:)';
        v(b) = mean(rec,'omitnan');
    end
end

function d = pick(nm,k)
    switch nm, case 'full', d=k; case 'readout', d=1; case 'null', d=k-1; end
end
function out = ternary(c,a,b), if c, out=a; else, out=b; end, end
