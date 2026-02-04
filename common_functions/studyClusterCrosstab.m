function OUT = studyClusterCrosstab(studyLabel, subjectClusterLabels, varargin)
% studyClusterCrosstab  Contingency / permutation test for study vs cluster membership
%
% OUT = studyClusterCrosstab(studyLabel, subjectClusterLabels)
% OUT = studyClusterCrosstab(..., 'nperm',1000, 'targetCluster', 1, 'verbose',true)
%
% Inputs:
%  - studyLabel:           N x 1 vector of integers (study index for each subject)
%  - subjectClusterLabels: N x 1 vector of integer cluster labels (1..K) for each subject
%
% Optional name/value:
%  'nperm'        - number of permutations for permutation p-value (default 1000)
%  'targetCluster'- integer cluster id to treat as "positive" when returning binary group (default 1)
%  'verbose'      - true/false (default true)
%
% Output (struct OUT):
%  - OUT.studyIDs           : unique study IDs (Ns x 1)
%  - OUT.contingencyCounts  : Ns x K matrix of counts (rows=studies, cols=clusters)
%  - OUT.contingencyProps   : Ns x K matrix of proportions per study (rows sum to 1)
%  - OUT.assignedGroup      : Ns x 1 assigned cluster per study (majority vote; ties -> NaN)
%  - OUT.binaryGroup        : Ns x 1 logical = (assignedGroup == targetCluster)
%  - OUT.chi2stat           : observed chi2 statistic
%  - OUT.chi2df             : degrees of freedom
%  - OUT.chi2p              : chi2 p-value (asymptotic)
%  - OUT.cramersV           : Cramer's V effect size
%  - OUT.perm_p             : permutation p-value for chi2 (NaN if nperm==0)
%  - OUT.perm_stats         : vector of permutation chi2 stats (length nperm)
%  - OUT.table_text         : small text summary for quick printing
%
% Example:
%  OUT = studyClusterCrosstab(studyLabel, stats.best_cluster_labels, 'nperm',2000);

% ---- parse inputs ----
p = inputParser;
addParameter(p,'nperm',1000,@(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(p,'targetCluster',1,@(x)isnumeric(x) && isscalar(x));
addParameter(p,'verbose',true,@islogical);
parse(p,varargin{:});
nperm = p.Results.nperm;
targetCluster = p.Results.targetCluster;
verbose = p.Results.verbose;

% ---- basic checks ----
if numel(studyLabel) ~= numel(subjectClusterLabels)
    error('studyLabel and subjectClusterLabels must have same length.');
end
studyLabel = studyLabel(:);
subjectClusterLabels = subjectClusterLabels(:);

% unique studies and clusters
[studyIDs, ~, studyIdx] = unique(studyLabel, 'stable'); % studyIdx maps obs -> row
clusterIDs = unique(subjectClusterLabels);
K = numel(clusterIDs);
Ns = numel(studyIDs);

% build contingency counts (Ns x K)
counts = zeros(Ns, K);
for k = 1:K
    mask = subjectClusterLabels == clusterIDs(k);
    c = accumarray(studyIdx(mask), 1, [Ns 1]);
    counts(:, k) = c;
end

% proportions (per-study)
props = counts ./ (sum(counts,2) + eps);  % Ns x K; small eps to avoid div-by-zero

% majority-vote assignment per study
assigned = nan(Ns,1);
for s = 1:Ns
    row = counts(s,:);
    if all(row==0)
        assigned(s) = NaN;
    else
        [mval, midx] = max(row);
        % check tie
        ties = find(row == mval);
        if numel(ties) > 1
            assigned(s) = NaN; % tie
        else
            assigned(s) = clusterIDs(midx);
        end
    end
end

% binary grouping
binaryGroup = false(Ns,1);
binaryGroup(~isnan(assigned)) = assigned(~isnan(assigned)) == targetCluster;

% ---- overall chi-square test (studies x clusters) ----
obs = counts;         % observed counts
Ntotal = sum(obs(:));
rowSums = sum(obs,2);
colSums = sum(obs,1);
expct = (rowSums * colSums) / Ntotal; % outer product
% handle zeros in expected: exclude cells with exp==0 from chi2 sum (they contribute 0)
validMask = expct > 0;
chi2stat = sum(((obs(validMask) - expct(validMask)).^2) ./ expct(validMask));
df = (Ns - 1) * (K - 1);
chi2p = 1 - chi2cdf(chi2stat, df);

% Cramer's V
cramersV = sqrt( chi2stat / (Ntotal * min(Ns-1, K-1)) );

% ---- permutation test (if requested) ----
permStats = NaN(nperm,1);
perm_p = NaN;
if nperm > 0
    if verbose, fprintf('Running %d permutations for chi2 null... ', nperm); tic; end
    rngstate = rng; % preserve RNG
    parfor_or_for = false;
    try
        % try parfor if Parallel Toolbox available
        pct = isempty(gcp('nocreate'));
        if ~pct
            parfor_or_for = true;
        else
            parfor_or_for = false;
        end
    catch
        parfor_or_for = false;
    end
    % Use simple loop (parfor disabled unless user has parallel pool)
    for ip = 1:nperm
        % permute cluster labels across subjects (keeps study sizes)
        permLabels = subjectClusterLabels(randperm(numel(subjectClusterLabels)));
        % recompute contingency
        pc = zeros(Ns,K);
        for k = 1:K
            mask = permLabels == clusterIDs(k);
            pc(:,k) = accumarray(studyIdx(mask), 1, [Ns 1]);
        end
        % chi2 stat for permuted
        pN = sum(pc(:));
        rS = sum(pc,2); cS = sum(pc,1);
        e = (rS * cS) / (pN + eps);
        vm = e > 0;
        permStats(ip) = sum(((pc(vm) - e(vm)).^2) ./ e(vm));
    end
    rng(rngstate); % restore RNG
    if verbose, fprintf('done (%.1fs)\n', toc); end
    perm_p = (sum(permStats >= chi2stat) + 1) / (nperm + 1); % add-one correction
end

% text summary
ttext = sprintf('Chi2=%0.3g, df=%d, p_asym=%0.4g, perm_p=%s, CramerV=%0.3g', ...
    chi2stat, df, chi2p, ...
    ternary(@(x) sprintf('%0.4g',x), @(x) x, perm_p, 'NaN'), cramersV);

% package outputs
OUT.studyIDs = studyIDs;
OUT.contingencyCounts = counts;
OUT.contingencyProps = props;
OUT.assignedGroup = assigned;
OUT.binaryGroup = binaryGroup;
OUT.chi2stat = chi2stat;
OUT.chi2df = df;
OUT.chi2p = chi2p;
OUT.cramersV = cramersV;
OUT.perm_p = perm_p;
OUT.perm_stats = permStats;
OUT.table_text = ttext;

% optional printing
if verbose
    fprintf('Contingency table: studies x clusters = %d x %d\n', Ns, K);
    disp('Counts (first few studies):'); disp(array2table(counts(1:min(10,Ns),:), 'VariableNames', cellfun(@(x) sprintf('C%d',x), num2cell(clusterIDs), 'Uni',false)));
    fprintf('%s\n', ttext);
end

end

% small helper: ternary utility
function out = ternary(formatF, valF, x, fallback)
% ternary chooses formatting function if x is numeric scalar; else fallback string
if ~isnan(x)
    out = formatF(x);
else
    out = fallback;
end
end
