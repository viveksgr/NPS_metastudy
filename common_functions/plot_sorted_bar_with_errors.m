function plot_sorted_bar_with_errors(data, column_names)
% data: matrix with size [max_subjects x 11], each column has NaNs for missing subjects
% column_names: cell array of strings (1x11), names of the columns

% Compute mean and standard error, ignoring NaNs
means = nanmean(data);
stderr = nanstd(data) ./ sqrt(sum(~isnan(data)));

% Sort in decreasing order
[sorted_means, idx] = sort(means, 'descend');
sorted_stderr = stderr(idx);
sorted_data = data(:, idx);
sorted_names = column_names(idx);

% Bar plot
% figure; 
figure('Position',[0.5 0.5 640 480])
hold on;
b = bar(1:length(sorted_means), sorted_means, 'FaceColor', [0.6 0.6 0.9]);

% Error bars
errorbar(1:length(sorted_means), sorted_means, sorted_stderr, ...
         'k', 'linestyle', 'none', 'LineWidth', 1.2);

% Overlay individual data points as black dots
for i = 1:length(sorted_means)
    y = sorted_data(:, i);
    y = y(~isnan(y));
    x = repmat(i, size(y));
    scatter(x, y, 10, 'k', 'filled', 'jitter','on', 'jitterAmount',0.1);
end

% Formatting
xticks(1:length(sorted_means));
safe_labels = strrep(sorted_names, '_', '\_');
xticklabels(safe_labels);
% xticklabels(sorted_names);
xtickangle(45);
ylabel('Score');
title('Sorted Group Scores with Individual Data');
box on; grid on;

end
