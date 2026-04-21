function fig = create_wm_fig(data_cell,data_cell2,data_cell3,Group_labels)

% Extract white and gray matter
fun_wm = @(X) extract_gray_white_csf(X, 'masks', ...
    {'gray_matter_mask.nii', 'canonical_white_matter_thrp5_ero1.nii', ...
    'canonical_ventricles_thrp5_ero1.nii'});
% nstud = vertcat(cellfun(@(x) size(x.dat,2),data_cell));
data_cell_wm = cellfun(fun_wm, data_cell, 'UniformOutput', false);
data_wm = vertcat(data_cell_wm{:});
data_cell_wm = cellfun(fun_wm, data_cell2, 'UniformOutput', false);
data_wm2 = vertcat(data_cell_wm{:});
data_cell_wm = cellfun(fun_wm, data_cell3, 'UniformOutput', false);
data_wm3 = vertcat(data_cell_wm{:});

nstud = vertcat(cellfun(@(x) size(x.dat,2),data_cell));
[tbl, X, groupVec, varNames] = createMixedDesign(Group_labels, nstud);
uniqueGroups = categories(categorical(Group_labels));
colorIdx = zeros(length(data_wm(:,1)),3);
cl_mat = lines(length(unique(Group_labels)));
for i = 1:length(data_wm(:,1))
    % find which canonical group this row belongs to (stable mapping)
    colorIdx(i,:) = cl_mat(find(strcmp(string(uniqueGroups), string(tbl{i,2}))),:);
end

% Plot
fig = figure('Position',[0.5 0.5 640 480]);
hold on
subplot(2,3,1)
plot_scatter_linear(data_wm(:,3),data_wm(:,1),colorIdx)
xlabel('CSF')
ylabel('GM')
% title('Raw Data')

subplot(2,3,4)
plot_scatter_linear(data_wm(:,2),data_wm(:,1),colorIdx)
xlabel('WM')
ylabel('GM')
% title('Raw Data')

subplot(2,3,2)
plot_scatter_linear(data_wm3(:,3),data_wm3(:,1),colorIdx)
xlabel('CSF')
ylabel('GM')
% title('renormed (add only)')

subplot(2,3,5)
plot_scatter_linear(data_wm3(:,2),data_wm3(:,1),colorIdx)
xlabel('WM')
ylabel('GM')
% title('renormed (add only)')

subplot(2,3,3)
plot_scatter_linear(data_wm2(:,3),data_wm2(:,1),colorIdx)
xlabel('CSF')
ylabel('GM')
% title('renormed (add+mul)')

subplot(2,3,6)
plot_scatter_linear(data_wm2(:,2),data_wm2(:,1),colorIdx)
xlabel('WM')
ylabel('GM')
% title('renormed (add+mul)')
