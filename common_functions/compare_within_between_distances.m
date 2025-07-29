function compare_within_between_distances(data)
% data: [subjects x conditions], NaNs allowed

[nSubj, nCond] = size(data);

within_diffs = [];
between_diffs = [];

for i = 1:nCond
    col_data = data(:, i);
    valid_idx = find(~isnan(col_data));
    vals = col_data(valid_idx);
    
    % All pairwise differences within this column
    for j = 1:length(vals)
        for k = j+1:length(vals)
            within_diffs(end+1) = abs(vals(j) - vals(k));
        end
    end
end

% Between-group diffs
for i = 1:nCond
    for j = i+1:nCond
        col1 = data(:, i);
        col2 = data(:, j);
        idx1 = find(~isnan(col1));
        idx2 = find(~isnan(col2));
        for a = 1:length(idx1)
            for b = 1:length(idx2)
                between_diffs(end+1) = abs(col1(idx1(a)) - col2(idx2(b)));
            end
        end
    end
end

% Stats and plot
fprintf('Within-mean = %.3f | Between-mean = %.3f\n', mean(within_diffs), mean(between_diffs));
[h,p] = ttest2(within_diffs, between_diffs);
fprintf('T-test: t=%.3f, p=%.3g\n', h, p);

figure;
boxplot([within_diffs'; between_diffs'], ...
        [ones(length(within_diffs),1); 2*ones(length(between_diffs),1)], ...
        'Labels', {'Within-group', 'Between-group'});
ylabel('Absolute Score Difference');
title('Comparison of Within vs. Between Subject Differences');
grid on;
end
