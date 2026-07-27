function R = run_mixture_withinsubj(varargin)
% RUN_MIXTURE_WITHINSUBJ  Within-subject validation of the mixture model on
% the crossover studies (5, 7, 9). The same subject appears in two arms on
% the same scanner, so person traits AND acquisition cancel.
%
% For each crossover study S:
%   - train the 8-class mixture model (whitened-space LDA, uniform priors)
%     and the readout w on ALL OTHER studies (S fully held out)
%   - project both arms of each subject -> per-subject 2 x 8 posterior profile
%   - A = lower-relief arm, B = higher-relief arm (so A->B increases relief)
%
% Metrics
%   (1) RESPONSIVENESS (paired, within-subject):
%       R_k = [ (P(A|A-state)+P(B|B-state)) - (P(A|B-state)+P(B|A-state)) ]/2
%       R_k>0  <=>  each state loads more on its own arm than the other.
%       Report mean R, signrank p, and within-subject ordering accuracy.
%   (2) SWITCH VECTOR decomposition:
%       mean switch dz = mean_k ( z_k^B - z_k^A ), split into on-w and null.
%       Sign-flip permutation tests whether the direction is consistent.
%       The null-space component = "different coordinates, not just more relief".
%   (3) TRAIT STABILITY: within-subject corr between the two arms' full
%       profiles (high => person carries a stable mixture across intervention).
%
% Usage:
%   R = run_mixture_withinsubj();                       % COORDS {C1,C2,C3}
%   R = run_mixture_withinsubj('COORDS',{'C1','C2','C3','NPS','SIIPS'});

p = inputParser;
addParameter(p,'MAT','datamat.mat');
addParameter(p,'COORDS',{'C1','C2','C3'});
addParameter(p,'N_PERM',5000);
addParameter(p,'RNG_SEED',7);
parse(p,varargin{:});
cfg = p.Results; rng(cfg.RNG_SEED);

GRP = {'Cognitive','Conditioning','Control','Mindfulness', ...
       'Placebo','Placebo_C','Remifentanil','Social'};
C = cfg.COORDS;

S = load(cfg.MAT); T = S.tbl2;
relief = -T.pain;
D = T{:,GRP}; [~,gi] = max(D,[],2);
g   = categorical(gi,1:numel(GRP),GRP);
sid = categorical(T.StudyID);
X   = T{:,C};

fprintf('=== WITHIN-SUBJECT MIXTURE VALIDATION | COORDS={%s} ===\n\n', strjoin(C,','));

% crossover studies = studies with >1 arm
studies = categories(removecats(sid)); cross = {};
for i = 1:numel(studies)
    if numel(categories(removecats(g(sid==studies{i})))) > 1, cross{end+1} = studies{i}; end %#ok<AGROW>
end

R = struct('study',{},'A',{},'B',{},'n',{},'meanR',{},'accR',{},'pR',{}, ...
           'onaxis',{},'nullmag',{},'reliefdiff',{},'p_on',{},'p_null',{}, ...
           'traitr',{},'profA',{},'profB',{},'names',{}, ...
           'diagv',{},'offv',{},'trr',{},'Rk',{},'pBB',{},'pBA',{},'pAA',{},'pAB',{}, ...
           'PA',{},'PB',{},'cA',{},'cB',{},'rA',{},'rB',{});

for c = 1:numel(cross)
    Sid = cross{c};
    te = (sid==Sid); tr = ~te;

    % train-only standardization + whitening + readout
    mu = mean(X(tr,:),1); sg = std(X(tr,:),0,1); sg(sg==0)=1;
    Xtr = (X(tr,:)-mu)./sg;  Xte = (X(te,:)-mu)./sg;
    Sig = resid_cov(Xtr, g(tr), sid(tr)); Wm = inv_sqrtm(Sig);
    Ztr = Xtr*Wm; Zte = Xte*Wm;
    what = fit_readout(Ztr, relief(tr), sid(tr), C);

    % mixture model: whitened-space LDA, uniform priors
    gtr = removecats(g(tr));
    mdl = fitcdiscr(Ztr, gtr, 'DiscrimType','pseudoLinear', 'Prior','uniform');
    [~, post] = predict(mdl, Zte);                 % nTe x nClass posteriors
    cls = string(mdl.ClassNames);

    % arms, ordered A=low relief, B=high relief
    gte = removecats(g(te)); rte = relief(te);
    arms = categories(gte);
    mr = [mean(rte(gte==arms{1})), mean(rte(gte==arms{2}))];
    if mr(1) <= mr(2), A = arms{1}; B = arms{2}; else, A = arms{2}; B = arms{1}; end

    idxA = find(gte==A); idxB = find(gte==B);
    n = min(numel(idxA), numel(idxB));
    idxA = idxA(1:n); idxB = idxB(1:n);             % index-aligned subjects

    cA = find(cls==A); cB = find(cls==B);
    postA = post(idxA,:); postB = post(idxB,:);
    pAA = postA(:,cA); pAB = postA(:,cB);
    pBB = postB(:,cB); pBA = postB(:,cA);

    % (1) responsiveness
    Rk  = ((pAA + pBB) - (pAB + pBA))/2;
    meanR = mean(Rk); accR = mean(Rk>0);
    pR  = signrank_safe(Rk);

    % (2) switch vector decomposition (B - A), whitened
    dz  = Zte(idxB,:) - Zte(idxA,:);
    mbar = mean(dz,1)';
    onax = what'*mbar; nullv = mbar - onax*what; nullmag = norm(nullv);
    reliefdiff = mean(rte(idxB)) - mean(rte(idxA));
    % sign-flip permutation: consistent direction?
    pon = 0; pnu = 0;
    for r = 1:cfg.N_PERM
        s = (rand(n,1)<0.5)*2-1;
        mb = mean(dz.*s,1)';
        o = what'*mb; nu = norm(mb - o*what);
        pon = pon + (abs(o) >= abs(onax));
        pnu = pnu + (nu >= nullmag);
    end
    p_on = (1+pon)/(cfg.N_PERM+1); p_null = (1+pnu)/(cfg.N_PERM+1);

    % (3) trait stability: within-subject profile correlation across arms
    tr_r = nan(n,1);
    for k = 1:n
        tr_r(k) = corr(postA(k,:)', postB(k,:)');
    end
    traitr = mean(tr_r,'omitnan');

    diagv = (pAA + pBB)/2;      % per-subject loading on the arm they are in
    offv  = (pAB + pBA)/2;      % per-subject loading on the other arm

    R(end+1) = struct('study',Sid,'A',A,'B',B,'n',n,'meanR',meanR,'accR',accR, ...
        'pR',pR,'onaxis',onax,'nullmag',nullmag,'reliefdiff',reliefdiff, ...
        'p_on',p_on,'p_null',p_null,'traitr',traitr, ...
        'profA',mean(postA,1),'profB',mean(postB,1),'names',{cls}, ...
        'diagv',diagv,'offv',offv,'trr',tr_r,'Rk',Rk, ...
        'pBB',pBB,'pBA',pBA,'pAA',pAA,'pAB',pAB, ...
        'PA',postA,'PB',postB,'cA',cA,'cB',cB, ...
        'rA',rte(idxA),'rB',rte(idxB)); %#ok<AGROW>

    fprintf('--- study %s | %s (A, low relief) <-> %s (B, high relief) | n=%d ---\n', Sid, A, B, n);
    fprintf('(1) responsiveness: meanR=%+.3f  ordering acc=%.2f  signrank p=%.4f\n', meanR, accR, pR);
    fprintf('    P(A|A)=%.3f P(A|B)=%.3f | P(B|B)=%.3f P(B|A)=%.3f\n', ...
        mean(pAA), mean(pBA), mean(pBB), mean(pAB));
    fprintf('(2) switch B-A: on-axis=%+.3f (p=%.4f) | null=%.3f (p=%.4f) | relief diff=%+.2f\n', ...
        onax, p_on, nullmag, p_null, reliefdiff);
    fprintf('(3) trait stability (within-subj profile corr) = %+.3f\n', traitr);
    % off-diagonal: what does arm B load on besides B?
    [~,ord] = sort(mean(postB,1),'descend');
    fprintf('    mean B-state profile top-3: ');
    for q=1:3, fprintf('%s=%.2f ', cls(ord(q)), mean(postB(:,ord(q)))); end
    fprintf('\n\n');
end

% aggregate responsiveness across studies (Stouffer on signrank not trivial; report simple)
fprintf('=== summary ===\n');
fprintf('%-8s %-14s %-14s %5s %7s %5s %6s %8s %8s\n', ...
    'study','A(low)','B(high)','n','meanR','acc','pR','null_p','trait_r');
for i = 1:numel(R)
    fprintf('%-8s %-14s %-14s %5d %+7.3f %5.2f %6.3f %8.4f %+8.3f\n', ...
        R(i).study, R(i).A, R(i).B, R(i).n, R(i).meanR, R(i).accR, R(i).pR, R(i).p_null, R(i).traitr);
end
end

% ================= helpers =================
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
function what = fit_readout(Z, y, sid, C)
    W = table(); W.relief = y; W.StudyID = removecats(categorical(sid));
    for j = 1:size(Z,2), W.(C{j}) = Z(:,j); end
    lme = fitlme(W, sprintf('relief ~ %s + (1|StudyID)', strjoin(C,' + ')));
    fe = lme.Coefficients; w = zeros(size(Z,2),1);
    for j = 1:size(Z,2), w(j) = fe.Estimate(strcmp(fe.Name,C{j})); end
    what = w/norm(w);
end
function p = signrank_safe(x)
    x = x(~isnan(x));
    if numel(x) < 2 || all(x==0), p = 1; return; end
    try, p = signrank(x); catch, [~,p] = ttest(x); end
end
