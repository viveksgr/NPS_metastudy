%% JKCorr stuff
% [sim_values, d, low_agreement, Nvox] = jackknife_similarity(data_cell_rn{1});
[sim_values, ~, d] = cellfun(@(x) jackknife_similarity(x),data_cell2,'UniformOutput',false);
% [sim_values_sad, ~, d1] = cellfun(@(x) jackknife_similarity(x,'similarity_metric','standardized_abs_deviation'),data_cell_rn2,'UniformOutput',false);
% [sim_values_msz, ~, d2] = cellfun(@(x) jackknife_similarity(x,'similarity_metric','mean_shift_z'),data_cell_rn2,'UniformOutput',false);
% [sim_values_ssz, ~, d3] = cellfun(@(x) jackknife_similarity(x,'similarity_metric','scale_shift_z'),data_cell_rn2,'UniformOutput',false);
% % 
% global_mean_sim = mean(vertcat(sim_values{:}));
% sim_values2 = cellfun(@(x) x-mean(x),sim_values,'UniformOutput',false);
% sim_values = cellfun(@(x) x-global_mean_sim,sim_values,'UniformOutput',false);

d_vec = ~vertcat(d{:});
v_vec = vertcat(sim_values2{:});
v_vec = v_vec-min(v_vec);
[tmap_iv3, pmap_iv, df,ts] = voxelwiseLM(data_cell_rn2, C, st_vec, 1,[],v_vec);


s_val_mat = cohensD;
s_val_mat.d = cellfun(@(x) mean(x),sim_values)';
s_val_mat.SE = cellfun(@(x) std(x)/sqrt(length(x)),sim_values)';
[r,p] = corrcoef(s_val_mat.d ,cohensD.d);
str = sprintf('Corr with NPS: %.3f. p: %.3f',r(2),p(2));
[figOut, ~] = plotCohensD_byGroup(s_val_mat, Group_labels, ...
    ...         % tells function to plot into this subplot
    'Colors', cl_mat, ...        % supply full color matrix matching canonical labels
    'Sorting', true, ...
    'ShowLegend', true,...
    'ShowXlabel', true ,...
    'GroupGap',1,...
    'FigureTitle',str);

sim_values_ = sim_values;
st=plotCellScatterFits(sim_values_, NPS_cell);
gcf
xlabel('JKCorr')
ylabel('NPS')
title(sprintf('m: %.2f; t: %.2f; p: %.2f',st.meanSlope,st.ttest.tstat, st.ttest.p))

%GLM
sim_vect = vertcat(sim_values_{:});
NPS_sc_reg = regress_out(NPS_sc, sim_vect);
%% Older code





%% Thresholds and counts

tmpl = image_vector('image_names', 'fmriprep20_template.nii');
tiv_rs = resample_space(tmap_iv, tmpl);
T = tinv(1 - p1, df);
plot_signed_threshold(tiv_rs, T)


[meanCountp, semCountp,rcp,rcpm] = cellfun(@(x) countThresholdedROIs(x, idx_pos),contrastCells,'UniformOutput',false);
[meanCountn, semCountn,rcn,rcnm] = cellfun(@(x) countThresholdedROIs(x, idx_neg),contrastCells,'UniformOutput',false);

s_val_mat = cohensD;
s_val_mat.d = vertcat(meanCountp{:});
s_val_mat.SE = vertcat(semCountp{:});
str = 'mean pos ROI';
[figOut, ordp] = plotCohensD_byGroup(s_val_mat, Group_labels, ...
    ...         % tells function to plot into this subplot
    'Colors', cl_mat, ...        % supply full color matrix matching canonical labels
    'Sorting', true, ...
    'ShowLegend', true,...
    'ShowXlabel', true ,...
    'GroupGap',0,...
    'FigureTitle',str);

s_val_mat = cohensD;
s_val_mat.d = vertcat(meanCountn{:});
s_val_mat.SE = vertcat(semCountn{:});
str = 'mean neg ROI';
[figOut, ~] = plotCohensD_byGroup(s_val_mat(ordp,:), Group_labels(ordp), ...
    ...         % tells function to plot into this subplot
    'Colors', cl_mat, ...        % supply full color matrix matching canonical labels
    'Sorting', false, ...
    'ShowLegend', true,...
    'ShowXlabel', true ,...
    'GroupGap',0,...
    'FigureTitle',str);

% Scatter 
uniqueGroups = categories(categorical(Group_labels));
colorIdx = zeros(nstud,3);
for i = 1:nstud
    % find which canonical group this row belongs to (stable mapping)
    colorIdx(i,:) = cl_mat(find(strcmp(string(uniqueGroups), string(Group_labels(i)))),:);
end
plot_scatter_linear(vertcat(meanCountp{:}),vertcat(meanCountn{:}),colorIdx)
plot(meanCountp',meanCountp','k')

tst = plotCellScatterFits(rcp, rcn);

%Study matrices
rcpm_m = cellfun(@(x) mean(x),rcnm,'UniformOutput',false);
rcpm_m = vertcat(rcpm_m{:});
% rcpm_m = rcpm_m-mean(rcpm_m,);

figure()
hold on
idxc = idx_neg;
lbls = cellfun(@(x) x(1:min(30,length(x))),atlas.label_descriptions(idxc),'UniformOutput',false);
imagesc(rcpm_m)
xticks(1:length(idxc))
xticklabels(lbls)
yticks(1:nstud)
yticklabels(ColNames)
axis tight
colorbar
hold on
ax = gca;
ax.XAxis.TickLabelInterpreter = 'none';
ax.YAxis.TickLabelInterpreter = 'none';
clim([-1 1])


%% LDA and decoding
addpath(genpath('C:\Work\Toolboxes\CanlabCore\custom_tools\roi_lda_decoding'));

results = canlab_roiwise_lda_decode(data_cell_rn2, Group_labels, ...
    'atlas', 'canlab2024_coarse_fmriprep20_2mm', ...
    'nFolds', 10, ...
    'leaveStudyOut', false, ...
    'nPermutations', 100, ...
    'useParallel', true, ...
    'parallelMode', 'permutations', ...
    'doplot', false, ...
    'saveModels', false);


% results = canlab_roiwise_lda_decode(data_cell_rn2, Group_labels, ...
%     'atlas', 'canlab2024_coarse_fmriprep20_2mm', ...
%     'nFolds', 2, ...
%     'leaveStudyOut', false, ...
%     'nPermutations', 100, ...
%     'randomSeed', 1, ...
%     'doplot', true);
 
% 
% results = canlab_roiwise_lda_decode(data_cell_rn2, Group_labels, ...
%     'atlas', 'canlab2024_coarse_fmriprep20_2mm', ...
%     'nFolds', 5, ...
%     'leaveStudyOut', true, ...
%     'doplot', false);

% tmap_iv2  =tmap_iv;

tmap_iv2 = results.p_value_map;
tmap_iv2.dat = 1- tmap_iv2.dat; 
% 
tmap_iv2 = results.overall_accuracy_map;
tmap_iv2.dat = tmap_iv2.dat/100; 
tmpl = image_vector('image_names','fmriprep20_template.nii');

% tmpl = image_vector('image_names','gray_matter_mask.nii');
tiv_rs = resample_space(tmap_iv2, tmpl);
plot_signed_threshold(tiv_rs, 0.975)


%%

% plot(meanCountp',0.066-meanCountp','k')

fun_demean = @(X) X-mean(X,'all'); 
fun = @(X) X(:,idx);
% [V_out, ~, ~] = regress_out_wm_csf(data_cell_rn{2}, nuisance_cell{2}, 'addIntercept',false);

contrastCells_red_ = cellfun(fun_demean, contrastCells, 'UniformOutput', false);
% contrastCells_red = cellfun(fun, contrastCells_red_, 'UniformOutput', false);
% plotInterventionROIHeat(contrastCells_red)

fun_avgsub = @(y) mean(y);
contrastCells_y = cellfun(fun_avgsub, contrastCells_red_, 'UniformOutput', false);
cm = vertcat(contrastCells_y{:});

opts.nComp = 10;
opts.doPCA = true;
opts.nPC   = 50;
opts.icaMethod = 'fastica';
[studyWeights, A_mix, S_maps, subj2study] = study_level_ICA( contrastCells_red_, opts);
% then cluster studyWeights (S x P)

% % M is nROI x nInterventions
% [fig, order, S] = spectralReorderAndPlot(cm', 'Labels', cohensD.Name, 'Affinity','corr', 'Title','ROI x Interventions (spectral reorder)');
% 
% % M is nStudy x nROI
% [fig, order, clusterIdx] = heatmap_sorted_rows(cm, 'RowLabels', cohensD.Name, 'ColLabels', atlas.labels(idx));

stats = clusterdata_permtest(cm, ...
    'k', 2:3, ...
    'distancemetric', 'correlation', ...
    'linkagemethod', 'average', ...
    'reducedims', false, ...
    'nperm', 100, ...
    'doplot', true, ...
    'verbose', true);

% Xticks
gca
x2 = str2num(xticklabels);
[~,argsort ] = sort(x2);
ax = gca;
ax.TickLabelInterpreter = 'none';
xticklabels(P_labels(argsort))

[cnum,argsort]=sort(stats.best_cluster_labels);
cm2 = cm(argsort,:);

% % Raw figure
imagesc(cm2(:,idx))
xticks(1:length(idx))
xticklabels(atlas.labels(idx))
yticks(1:size(cm,1))
% yticklabels(c2(argsort))
xline(sum(idx_pos)+0.5)
yline(sum(cnum==1)+0.5)
yline(sum(cnum<3)+0.5)
hold on
ax = gca;
ax.XAxis.TickLabelInterpreter = 'none';
ax.YAxis.TickLabelInterpreter = 'none';

% Cluster

%% GLM
C = stats.best_cluster_labels;
C(or(C==1,C==2))=0;
C(C>1)=1;
data_cell_rn2(del_cell)=[];
st_vec(del_cell)=[];
[tmap_iv, pmap_iv, df] = voxelwiseLM(data_cell_rn2, C, st_vec,2);


%% subject-wise c
[meanr, stdr] = cellfun(@(x) mean_subject_corr(x),data_cell_rn,'UniformOutput',true);

[mean_r_crs, se_r_crs, mean_z, se_z, nPairs] = mean_crossstudy_subject_corr(data_cell_rn);
rem_idx = 1:31;
cohensD_subjerr = cohensD;
cohensD_subjerr.d = meanr(rem_idx)';%'-mean_r_crs(rem_idx);
cohensD_subjerr.SE = stdr(rem_idx)';
[figOut, order] = plotCohensD_byGroup(cohensD_subjerr, P_labels, ...
     ...         % tells function to plot into this subplot
        'Colors', cl_mat, ...        % supply full color matrix matching canonical labels
        'Sorting', true, ...
        'ShowLegend', true,...
        'ShowXlabel', true ,...
        'GroupGap',2,...
        'FigureTitle',labels{1}); 

inc_idx = true(1,31);
% inc_idx(13) = false;
figure()
hold on
subplot(1,3,1)
plot_scatter_linear(cohensD_subjerr.d(inc_idx) ,P_pain(inc_idx),colorIdx(inc_idx,:))
xlabel('Cross subject correlation')
ylabel('pain scores')
subplot(1,3,2)
plot_scatter_linear(cohensD_subjerr.d(inc_idx) ,P_SIIPS(inc_idx),colorIdx(inc_idx,:))
xlabel('Cross subject correlation')
ylabel('SIIPS')
subplot(1,3,3)
plot_scatter_linear(cohensD_subjerr.d(inc_idx) ,P_NPS(inc_idx),colorIdx(inc_idx,:))
xlabel('Cross subject correlation')
ylabel('NPS')

%% subject-wise c



%% Common factors +ve -ve
% Apply parcellation
atlasname = 'canlab2024_coarse_fmriprep20_2mm';
atlas = load_atlas(atlasname);
f_correct = cell(1,31);
[contrastCells] = apply_parc_contrasiobj(data_cell2,atlas,f_correct);  

% contrastCells: 1 x Nc, each cell Ns_c x Nparcels
Nc = numel(contrastCells);
n_sub = cellfun(@(x) size(x,1), contrastCells);

% make big data matrix: rows = subjects, cols = parcels
data = vertcat(contrastCells{:});   % (sum Ns_c) x Nparcels

% create study label for each subject (useful for post-hoc summaries)
studyLabel = repelem(1:Nc, n_sub)'; % column, same length as rows in data

k_range = 2:8;
nperm = 500;  % raise to 1k+ for final tests

[cleanData, kept, removed, stats] = remove_low_variance_features(data, ...
    'RelSTDThreshold', 1e-3, 'AbsSTDThreshold', 1e-8, 'Names', atlas.label_descriptions);

stats = clusterdata_permtest(cleanData, ...
    'k', k_range, ...
    'distancemetric', 'correlation', ...
    'linkagemethod', 'average', ...
    'reducedims', true, ...
    'nperm', nperm, ...
    'doplot', true, ...
    'verbose', true);


%% Crosstab analysis
[H,~,leafOrder] = dendrogram(stats.linkage_tree, 0); % 
% how many leaves per study cluster in a given cluster partition?
bestLabels = stats.best_cluster_labels; % length = number of rows (subjects)

% your inputs:
% studyLabel: 826×1 integer vector mapping subject->study
% stats.best_cluster_labels: 826×1 integer vector of cluster labels (1 or 2 here)

OUT = studyClusterCrosstab(studyLabel, stats.best_cluster_labels, 'nperm', 2000, 'targetCluster', 1);
OUT2 = reorderContingencyRows(OUT, 'K', 2, 'Distance','correlation', 'Linkage','average', 'Plot', true);

% quick inspect:
imagesc(OUT2.contingencyProps_re)         % rows = studies, cols = clusters
xlabel('Cluster'); ylabel('Study (index)');
colorbar; title('Per-study proportion of subjects in each cluster');

% study-level binary grouping:
studyBinary = OUT.binaryGroup;    % length(unique(studyLabel))

%% Study maps
% % map subjects to clusters, but in dendrogram order:
% clustersInLeafOrder = bestLabels(leafOrder);
% tab = table(orderedStudyPerLeaf(:), clustersInLeafOrder(:), leafGroup(:));
% % a quick cross-tab
% crosstab(orderedStudyPerLeaf, clustersInLeafOrder)

% Clusterwise maps
[tmap_iv, pmap_iv, df] = voxelwiseLM(data_cell_rn(studyBinary), C(studyBinary), st_vec(studyBinary));


%% Raw correlation NPS and pain ratings
cohensD_NPS = cohensD_table_wrapper(Sgn_dat.NPS);
cohensD_SIIPS = cohensD_table_wrapper(Sgn_dat.SIIPS);
cohensD_rawpain = cohensD_table_wrapper(BigDat_n);
Group_labels = ColNames_mat.group_var;
rem_idx = true(1,length(Group_labels));
rem_idx([4 10 29])=false;

P_NPS = cohensD_NPS.d(rem_idx);
P_SIIPS = cohensD_SIIPS.d(rem_idx);
P_pain = cohensD_rawpain.d(rem_idx);
P_labels = Group_labels(rem_idx);
% P_labels = Group_labels2;
cl_mat = lines(length(unique(P_labels)));

uniqueGroups = categories(categorical(P_labels));
colorIdx = zeros(sum(rem_idx),3);
for i = 1:sum(rem_idx)
    % find which canonical group this row belongs to (stable mapping)
    colorIdx(i,:) = cl_mat(find(strcmp(string(uniqueGroups), string(P_labels(i)))),:);
end
plot_scatter_linear(P_NPS,P_pain,colorIdx)
scatter(P_NPS,P_SIIPS,36,colorIdx,'filled', 'MarkerEdgeColor', 'k')
% plot_scatter_linear(P_NPS,P_SIIPS,colorIdx)

% [fog, st] = plotResidualsByMedianGroup(P_NPS, P_pain); 
%% Atlas

% Apply parcellation
atlasname = 'canlab2024_coarse_fmriprep20_2mm';
atlas = load_atlas(atlasname);
f_correct = cell(1,30);
[contrastCells] = apply_parc_contrasiobj(data_cell_rn2,atlas,f_correct);  
contrastCells(del_cell)=[];

c_map = extract_roi_averages(tmap_iv, atlas);
roi_means = zeros(1,length(atlas.labels));
for nn=1:length(atlas.labels); roi_means(nn)=median(c_map(nn).dat); end

idx = 1:length(atlas.labels);
idx_pos = roi_means>ts(1);
idx_neg = roi_means<-ts(2);
idx = [idx(idx_pos) idx(idx_neg)];

fun_demean = @(X) X-mean(X,'all'); 
fun = @(X) X(:,idx);
% [V_out, ~, ~] = regress_out_wm_csf(data_cell_rn{2}, nuisance_cell{2}, 'addIntercept',false);

contrastCells_red_ = cellfun(fun_demean, contrastCells, 'UniformOutput', false);
% contrastCells_red = cellfun(fun, contrastCells_red_, 'UniformOutput', false);
% plotInterventionROIHeat(contrastCells_red)

fun_avgsub = @(y) mean(y);
contrastCells_y = cellfun(fun_avgsub, contrastCells_red_, 'UniformOutput', false);
cm = vertcat(contrastCells_y{:});

opts.nComp = 10;
opts.doPCA = true;
opts.nPC   = 50;
opts.icaMethod = 'fastica';
[studyWeights, A_mix, S_maps, subj2study] = study_level_ICA( contrastCells_red_, opts);
% then cluster studyWeights (S x P)

% % M is nROI x nInterventions
% [fig, order, S] = spectralReorderAndPlot(cm', 'Labels', cohensD.Name, 'Affinity','corr', 'Title','ROI x Interventions (spectral reorder)');
% 
% % M is nStudy x nROI
% [fig, order, clusterIdx] = heatmap_sorted_rows(cm, 'RowLabels', cohensD.Name, 'ColLabels', atlas.labels(idx));

stats = clusterdata_permtest(cm, ...
    'k', 2:3, ...
    'distancemetric', 'correlation', ...
    'linkagemethod', 'average', ...
    'reducedims', false, ...
    'nperm', 100, ...
    'doplot', true, ...
    'verbose', true);

% Xticks
gca
x2 = str2num(xticklabels);
[~,argsort ] = sort(x2);
ax = gca;
ax.TickLabelInterpreter = 'none';
xticklabels(P_labels(argsort))

[cnum,argsort]=sort(stats.best_cluster_labels);
cm2 = cm(argsort,:);

% % Raw figure
imagesc(cm2(:,idx))
xticks(1:length(idx))
xticklabels(atlas.labels(idx))
yticks(1:size(cm,1))
% yticklabels(c2(argsort))
xline(sum(idx_pos)+0.5)
yline(sum(cnum==1)+0.5)
yline(sum(cnum<3)+0.5)
hold on
ax = gca;
ax.XAxis.TickLabelInterpreter = 'none';
ax.YAxis.TickLabelInterpreter = 'none';

% Cluster

%% GLM
C = stats.best_cluster_labels;
C(or(C==1,C==2))=0;
C(C>1)=1;
data_cell_rn2(del_cell)=[];
st_vec(del_cell)=[];
[tmap_iv, pmap_iv, df] = voxelwiseLM(data_cell_rn2, C, st_vec,2);


%% subject-wise c
[meanr, stdr] = cellfun(@(x) mean_subject_corr(x),data_cell_rn,'UniformOutput',true);

[mean_r_crs, se_r_crs, mean_z, se_z, nPairs] = mean_crossstudy_subject_corr(data_cell_rn);
rem_idx = 1:31;
cohensD_subjerr = cohensD;
cohensD_subjerr.d = meanr(rem_idx)';%'-mean_r_crs(rem_idx);
cohensD_subjerr.SE = stdr(rem_idx)';
[figOut, order] = plotCohensD_byGroup(cohensD_subjerr, P_labels, ...
     ...         % tells function to plot into this subplot
        'Colors', cl_mat, ...        % supply full color matrix matching canonical labels
        'Sorting', true, ...
        'ShowLegend', true,...
        'ShowXlabel', true ,...
        'GroupGap',2,...
        'FigureTitle',labels{1}); 

inc_idx = true(1,31);
% inc_idx(13) = false;
figure()
hold on
subplot(1,3,1)
plot_scatter_linear(cohensD_subjerr.d(inc_idx) ,P_pain(inc_idx),colorIdx(inc_idx,:))
xlabel('Cross subject correlation')
ylabel('pain scores')
subplot(1,3,2)
plot_scatter_linear(cohensD_subjerr.d(inc_idx) ,P_SIIPS(inc_idx),colorIdx(inc_idx,:))
xlabel('Cross subject correlation')
ylabel('SIIPS')
subplot(1,3,3)
plot_scatter_linear(cohensD_subjerr.d(inc_idx) ,P_NPS(inc_idx),colorIdx(inc_idx,:))
xlabel('Cross subject correlation')
ylabel('NPS')

%% subject-wise c

[tbl, X, groupVec, varNames] = createMixedDesign(Group_labels, cohensD_subjerr.n, cohensD_subjerr.d, {'mean_corr'});

V = Sgn_dat.NPS(:,rem_idx);
V = V{:,:};
V = V(:);
V(isnan(V))=[];

tbl.Y = V;
lme = fitlme(tbl, 'Y ~ 1 + StudyLabel + mean_corr + (1|StudyID)');
disp(lme);


V = BigDat_n(:,rem_idx);
V = V(:);
V(isnan(V))=[];
tbl.Y = V;
lme = fitlme(tbl, 'Y ~ 1 + StudyLabel + mean_corr + (1|StudyID)');
disp(lme);