function crosstab_cluster_overlap(labels1, labels2, varargin)
% crosstab_cluster_overlap  Crosstab + % overlap barplot for two clustering vectors
%
% crosstab_cluster_overlap(labels1, labels2)
% crosstab_cluster_overlap(labels1, labels2, 'name1','ClustA', 'name2','ClustB')
%
% Inputs:
%   labels1, labels2 : N x 1 numeric vectors of cluster IDs (same length)
%
% Optional name/value:
%   'name1'   - label for first clustering (default: 'Clustering 1')
%   'name2'   - label for second clustering (default: 'Clustering 2')

p = inputParser;
addParameter(p, 'name1', 'Clustering 1', @ischar);
addParameter(p, 'name2', 'Clustering 2', @ischar);
parse(p, varargin{:});
name1 = p.Results.name1;
name2 = p.Results.name2;

labels1 = labels1(:);
labels2 = labels2(:);
assert(numel(labels1) == numel(labels2), 'labels1 and labels2 must be the same length.');

ids1 = unique(labels1);
ids2 = unique(labels2);
K1 = numel(ids1);
K2 = numel(ids2);

% Build contingency matrix (K1 x K2)
counts = zeros(K1, K2);
for i = 1:K1
    for j = 1:K2
        counts(i, j) = sum(labels1 == ids1(i) & labels2 == ids2(j));
    end
end

% Row-normalised proportions: for each cluster in labels1, % in each labels2 bin
props = counts ./ (sum(counts, 2) + eps);  % K1 x K2

% Chi-square test
[~, chi2p, chi2stat] = chi2cont(counts);

% Plot
figure('Color','w');
bar(props * 100, 'stacked');

tickLabels = arrayfun(@(x) sprintf('%d', x), ids1, 'UniformOutput', false);
xticks(1:K1);
xticklabels(tickLabels);
xlabel(name1);
ylabel('% of cluster');

legendLabels = arrayfun(@(x) sprintf('%s %d', name2, x), ids2, 'UniformOutput', false);
legend(legendLabels, 'Location', 'eastoutside');

title(sprintf('%s vs %s   \\chi^2=%.2f, p=%.3g', name1, name2, chi2stat, chi2p));

end


function [tbl, chi2stat, chi2p] = chi2cont(obs)
% Simple chi-square test of independence for a contingency matrix
Ntotal = sum(obs(:));
rowSums = sum(obs, 2);
colSums = sum(obs, 1);
expct = (rowSums * colSums) / Ntotal;
valid = expct > 0;
chi2stat = sum(((obs(valid) - expct(valid)).^2) ./ expct(valid));
df = (size(obs,1)-1) * (size(obs,2)-1);
chi2p = 1 - chi2cdf(chi2stat, df);
tbl = array2table(obs);
end
