% % gather all subject contrast files (example pattern)
% clear all
% confiles = dir(fullfile(pwd,'s_r_*.nii.gz')); % adjust
% % build fullpaths cell
% confiles = fullfile({confiles.folder}, {confiles.name});
% 
% % create fmri_data object (uses CanlabCore fmri_data constructor)
% DATA_OBJ_CON{1} = fmri_data(confiles);    % will apply default brain mask (or specify mask)
% save('contrast_data_objects.mat')

%% Simple plots
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
%%
QC metrics

% % c = qc_metrics_secon_level(BigDat_f{1}); 
% [outMat, fieldNames] = collateContrastScalars(BigDat_f,@(c) qc_metrics_second_level(c));
% figure('Position',[0.5 0.5 640 480])
% imagesc(outMat)
% [m,n] = size(outMat);
% xticks(1:1:n)
% xticklabels(ColNames_edit)
% xtickangle(90)
% yticks(1:1:m)
% yticklabels(fieldNames)
% ax = gca;
% ax.XAxis.TickLabelInterpreter = 'none';
% ax.YAxis.TickLabelInterpreter = 'none';
% 

% 
% 
% k = 0;
% figure('Position',[0.5 0.5 1280 720])
% hold on
% for ii = 17:31
%     k = k+1;
%     subplot(4,4,k)
%     plotGlobalMeans(data_cell_rn{ii})
%     ylim([-0.5 0.5])
%     title(ColNames_edit_full{ii},'Interpreter','none')
% end

%% Run fixed effects
plot_study_kmeans(dsc, 2, C )

computeContrastSimilarity(dsc,ColNames_edit);
[c2,argsort] = sort(C);
dsc2 = dsc(argsort);
ColNames2 = ColNames_edit(argsort);

computeContrastSimilarity(dsc2 ,ColNames2 );
computeContrastSimilarity_rsa(dsc2,ColNames2);

% Apply parcellation
atlasname = 'canlab2023_coarse_fmriprep20_1mm';
atlas = load_atlas(atlasname);
f_correct = cell(1,31);
[parc_data_multi] = apply_parc_contrasiobj(data_cell_rn,atlas,f_correct);  

parc_data_multi = parc_data_multi(argsort);
simMat = computeDistanceSimilarity_euc(parc_data_multi,ColNames2);

plotMDSFromSimilarity(simMat, c2, ColNames2)