addpath(genpath('C:\Work\Toolboxes'))

% Make simple barplot
ms_ratings_all_bc = ms_ratings_all;

% change pain_scores for positive effect`
sign_vector = [1 1 1 -1 1 1 -1 1 -1 -1 -1];

% pain_ratings:
ms_ratings_vec = sign_vector .*ms_ratings_all./std(ms_ratings_all,'omitmissing');

plot_sorted_bar_with_errors((ms_ratings_vec),study_names_correct)
compare_within_between_distances(ms_ratings_vec)
compare_avg_distances_per_column((ms_ratings_vec))

% Apply parcellation
atlasname = 'canlab2023_coarse_fmriprep20_1mm';
atlas = load_atlas(atlasname);
f_correct = cell(1,11);
f_correct{6} = [451   455   470   480   499   505   506];
f_correct{9} = [480];

% % Change sign
for ss = 1:nStudies
    ms_brain{ss}.dat = sign_vector(ss).* ms_brain{ss}.dat;
end
[parc_data_multi] = apply_parc_contrasiobj(ms_brain,atlas,f_correct);

% study_means = average_subjects(parc_data_neg);
[study_means] = average_subjects_t(parc_data_multi);

plot_top_roi_heatmap(study_means,atlas.labels, study_names, 1.65, 100)
plot_unique_roi_heatmap(  study_means, atlas.labels,  study_names)
[mean_data, var_meta] = average_subjects_tcorrected(parc_data_multi);

mean_ct_pos = sum(study_means>1.98,1);
mean_ct_neg = sum(study_means<-1.98,1);
create_nifty_patterns(atlas,-mean_data ,'mean_data_neg.nii')


group_idx = [1 1 1 1 2 1 1 1 2 2 2];
plot_study_kmeans(parc_data_multi, 2, group_idx)

%% Reintegrate data

studydir = 'C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\NPS_metastudy2.0\Data\subjectlevel\included_studies';
[BigDat, ColNames] = NPSMS_harvest_canlab_cols(studydir, ccode, 40);
BigDat = l2normalize_columns(BigDat,'L2');
BigDat = l2normalize_columns(BigDat,'InvVar');

% Switch signs
BigDat_n = BigDat;
idx = [3 4 7 13 14 15];
% ColNames_edit = ColNames; % Altered Manually
idx_vec = ones(1,size(BigDat,2));
idx_vec(idx) = -1;
idx_mat = repmat(idx_vec,[size(BigDat,1),1]);
BigDat_n = BigDat.*idx_mat;

plot_sorted_bar_with_errors(BigDat_n, ColNaphmes_edit)


%% Imaging Data
%  desc = descriptives(DATA_OBJ_CON{1, 1}  , ['noverbose', 'plotcoverage']);
%  stats = clusterdata_permtest(DATA_OBJ_CON{1, 1}.dat, 'k', [2:20], 'reducedims', true, 'ndims', 25);
%  % desc = descriptives(DATA_OBJ_CON{1, 1}  , ['noverbose', 'plotcoverage']);
%  % qc_metrics_second_level(DATA_OBJ_CON{1, 1});
% % help robfit_parcelwise

studydir = 'C:\Work\NPSSR\Data\subjectlevel';
[BigDat_f] = harvest_canlab_funccols(studydir, idx);

for ii = 1:16
    BigDat_f{ii}=  rescale(BigDat_f{ii}, 'l2norm_images');
end

% c = qc_metrics_second_level(BigDat_f{1});
[outMat, fieldNames] = collateContrastScalars(BigDat_f,@(c) qc_metrics_second_level(c));

figure('Position',[0.5 0.5 640 480])
imagesc(outMat)
[m,n] = size(outMat);
xticks(1:1:n)
xticklabels(ColNames)
xtickangle(90)
yticks(1:1:m)
yticklabels(fieldNames)
ax = gca;
ax.XAxis.TickLabelInterpreter = 'none';
ax.YAxis.TickLabelInterpreter = 'none';
k = 0;
figure('Position',[0.5 0.5 1280 720])
hold on
for ii = 1:n
    k = k+1;
    subplot(4,4,k)
    plotGlobalMeans(BigDat_f{ii})
    title(ColNames{ii},'Interpreter','none')
end

% Normalize to same space:
BigDat_n = alignFmriDataToReference(BigDat_f, [2]);



%% Run fixed effects