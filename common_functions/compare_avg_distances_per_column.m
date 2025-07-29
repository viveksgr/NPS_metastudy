function compare_avg_distances_per_column(data)
% data: [subjects x 11 columns], NaNs allowed

nCond = size(data, 2);
within_means = zeros(1, nCond);
between_means = zeros(1, nCond);

for i = 1:nCond
    % Within-group distances
    col_data = data(:, i);
    vals = col_data(~isnan(col_data));
    dists = [];
    for j = 1:length(vals)
        for k = j+1:length(vals)
            dists(end+1) = abs(vals(j) - vals(k));
        end
    end
    within_means(i) = mean(dists, 'omitnan');

    % Between-group distances
    between_dists = [];
    for j = 1:nCond
        if j == i, continue; end
        col_other = data(:, j);
        vals_j = col_other(~isnan(col_other));
        for a = 1:length(vals)
            for b = 1:length(vals_j)
                between_dists(end+1) = abs(vals(a) - vals_j(b));
            end
        end
    end
    between_means(i) = mean(between_dists, 'omitnan');
end

% T-test
fprintf('Mean within = %.3f | Mean between = %.3f\n', mean(within_means), mean(between_means));
[h, p, ci, stats] = ttest(within_means, between_means);
fprintf('Paired t-test: t=%.3f, p=%.4f\n', stats.tstat, p);

% Plot
figure;
boxplot([within_means', between_means'], ...
        [ones(nCond,1); 2*ones(nCond,1)], ...
        'Labels', {'Within', 'Between'});
ylabel('Mean Absolute Difference per Column');
title('Within vs Between Group Distances (Per Column)');
grid on; box on;
end
