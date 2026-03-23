addpath(genpath('C:\Work\Toolboxes'))
addpath('C:\Work\NPS_metastudy\common_functions')

%% Behavioral Datasets
% Non Placebo studies
ColNames_mat = load('C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\Postprocessing\labels_update.mat');
ColNames = ColNames_mat.ColNames;
Group_labels = ColNames_mat.Group_labels;

ccode = {[1],[3],[4],[21],[10,11,12],[15:17],[9],[19],[1],[3],[1]}';
studydir = 'C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\included_studies';
[BigDat, ~] = NPSMS_harvest_canlab_cols(studydir, ccode, 50);

% Placebo studies
studydir2 = 'C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\zunhammer2018_studies';
ccode2 = num2cell(ones(1,16));
[BigDat2, ~] = NPSMS_harvest_canlab_cols(studydir2, ccode2, 50);
BigDat = cat(2,BigDat, BigDat2);
lastRow = find(any(~isnan(BigDat),2), 1, 'last');
if ~isempty(lastRow); BigDat = BigDat(1:lastRow, :);end
% BigDat = l2normalize_columns(BigDat,'L2');
BigDat = l2normalize_columns(BigDat,'InvStd');

% Switch signs
BigDat_n = BigDat;
idx_vec = sign(mean(BigDat_n,'omitmissing'));
idx_mat = repmat(idx_vec,[size(BigDat,1),1]);
BigDat_n = -BigDat.*idx_mat;
plot_sorted_bar_with_errors(BigDat_n, ColNames)

%% Imaging Data
%  desc = descriptives(DATA_OBJ_CON{1, 1}  , ['noverbose', 'plotcoverage']);
%  stats = clusterdata_permtest(DATA_OBJ_CON{1, 1}.dat, 'k', [2:20], 'reducedims', true, 'ndims', 25);
%  % desc = descriptives(DATA_OBJ_CON{1, 1}  , ['noverbose', 'plotcoverage']);
%  % qc_metrics_second_level(DATA_OBJ_CON{1, 1});
%  % help robfit_parcelwise
idx_vec = [1 -1 1 -1 -1 1 -1 1 1 1 1 -1 -1 -1 1 ones(1,16)]';

st_vec_ind = cellfun(@(x) numel(x),ccode);
st_vec = [];
kk = 0;
for ii=1:length(st_vec_ind) 
    for jj = 1:st_vec_ind(ii)
        kk = kk+1;
        st_vec(kk)=ii;
    end
end

st_vec = cat(2,st_vec,12:27);
% Dummy contrast
C = (([1 1 1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 -1 -1 -1 1 1 -1 -1 -1 -1 -1 1 1 -1 1 1 -1]+1)/2)+1;

Idx = {[1], [1], [1], [2], [4:6],[1:3],[1],[1], [1], [1], [2]}';
[BigDat_f] = harvest_canlab_funccols(studydir, Idx);
ncontrast = numel(BigDat_f);
for ii = 1:ncontrast
    % BigDat_f{ii}=  rescale(BigDat_f{ii}, 'csf_mean_var');
    BigDat_f{ii}=  rescale(BigDat_f{ii}, 'l2norm_images');
    % BigDat_f{ii}=  rescale(BigDat_f{ii}, 'prctileimages');
    BigDat_f{ii}.dat  = -BigDat_f{ii}.dat.*idx_vec(ii); % Change signs of the contrasts so that all contrasts are for analgesia
end

% Zunhammer studies
Idx2 = num2cell(ones(16,1));
[BigDat_f2] = harvest_canlab_funccols(studydir2, Idx2);
ncontrast = numel(BigDat_f2);
for ii = 1:ncontrast
    % BigDat_f2{ii}=  rescale(BigDat_f2{ii}, 'csf_mean_var');
    BigDat_f2{ii}=  rescale(BigDat_f2{ii}, 'l2norm_images');
     % BigDat_f2{ii}=  rescale(BigDat_f2{ii}, 'prctileimages');
    BigDat_f2{ii}.dat  = BigDat_f2{ii}.dat; % Change signs of the contrasts so that all contrasts are for analgesia
end
BigDat_f = cat(2,BigDat_f, BigDat_f2);
% Normalize
data_cell = alignFmriDataToReference(BigDat_f, [2]);

%% Preprocessing/Cleanup

% Harmonize and negate
data_cell_rn = cellfun(@(x) harmonize_zero_preserve(x),data_cell,'UniformOutput',false);

% CSF
nuisance_cell = cellfun(@(x) extract_wm_csf_comps(x,false), data_cell_rn,'UniformOutput',false);
fun = @(V,X) regress_out_wm_csf(V, X, 'addIntercept',false,'verbose', false);
data_cell_rn2 = cellfun(fun, data_cell_rn, nuisance_cell, 'UniformOutput', false);
% CSF
% nuisance_cell = cellfun(@(x) extract_wm_csf_comps(x,false), data_cell,'UniformOutput',false);
% data_cell_2 = cellfun(fun, data_cell, nuisance_cell, 'UniformOutput', false);
[tmap_iv, pmap_iv, df,ts] = voxelwiseLM(data_cell_rn2, C, st_vec,1);

% Normalization
k = 0;
figure('Position',[0.5 0.5 1280 720])
hold on
for ii = 1:30
    k = k+1;
    subplot(5,6,k)
    plotGlobalMeans(data_cell_rn{ii})
    ylim([-0.5 0.5])
    title(ColNames{ii},'Interpreter','none')
end

% Remove wm
% [dx, dx2] = data_cell{1}.normalize_gm_by_wm_csf;

fun = @(X) normalize_gm_by_wm_csf(X);
data_cell_rn2 = cellfun(fun, data_cell_rn, 'UniformOutput', false);
data_cell2 = cellfun(fun, data_cell, 'UniformOutput', false);
data_cell2_rn = cellfun(@(x) harmonize_zero_preserve(x),data_cell2,'UniformOutput',false);


% WM mapping
fun_wm = @(X) extract_gray_white_csf(X, 'masks', ...
    {'gray_matter_mask.nii', 'canonical_white_matter_thrp5_ero1.nii', ...
    'canonical_ventricles_thrp5_ero1.nii'});
nstud = vertcat(cellfun(@(x) size(x.dat,2),data_cell_rn));
data_cell_wm = cellfun(fun_wm, data_cell_rn2, 'UniformOutput', false);
data_wm = vertcat(data_cell_wm{:});


[tbl, X, groupVec, varNames] = createMixedDesign(Group_labels, nstud);
uniqueGroups = categories(categorical(Group_labels));
colorIdx = zeros(length(data_wm(:,1)),3);
cl_mat = lines(length(unique(Group_labels)));
for i = 1:length(data_wm(:,1))
    % find which canonical group this row belongs to (stable mapping)
    colorIdx(i,:) = cl_mat(find(strcmp(string(uniqueGroups), string(tbl{i,2}))),:);
end
figure('Position',[0.5 0.5 640 480])
plot_scatter_linear(data_wm(:,3),data_wm(:,1),colorIdx)
xlabel('CSF')
ylabel('GM')
% legend(uniqueGroups)


% [outMat1, fieldNames] = collateContrastScalars(data_cell,@(c) qc_metrics_second_level(c));
% [outMat2, fieldNames] = collateContrastScalars(data_cell_rn_csf,@(c) qc_metrics_second_level(c));
% figure('Position',[0.5 0.5 640 480])
% subplot(2,1,1)
% imagesc(outMat1)
% [m,n] = size(outMat1);
% xticks([])
% % xticklabels(ColNames)
% xtickangle(90)
% yticks(1:1:m)
% yticklabels(fieldNames)
% ax = gca;
% ax.XAxis.TickLabelInterpreter = 'none';
% ax.YAxis.TickLabelInterpreter = 'none';
% subplot(2,1,2)
% imagesc(outMat2)
% [m,n] = size(outMat2);
% xticks(1:1:n)
% xticklabels(ColNames)
% xtickangle(90)
% yticks(1:1:m)
% yticklabels(fieldNames)
% ax = gca;
% ax.XAxis.TickLabelInterpreter = 'none';
% ax.YAxis.TickLabelInterpreter = 'none';

%% NPS masks
ord = load('C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\Postprocessing\ratings_sort.mat');
ord = ord.ord;
Sgn_dat = apply_all_signatures(data_cell2,'conditionnames',ColNames);
% Sgn_dat2 = apply_nps(data_cell,'conditionnames',ColNames);

cl_mat = lines(length(unique(Group_labels)));
% labels = {'NPS','NPSpos','NPSneg','SIIPS','PINES','Rejection','VPS','FM_Multisens','FM_pain','Empathic_Care'};
labels = {'NPS'};
ax = figure('Color','w','Position',[100 100 2400 1200]);
hold on
eval(sprintf('cohensD = cohensD_table_wrapper(Sgn_dat.%s);',labels{1}))

% cohensD = cohensD(ord,:);
% plotCohensD_byGroup(cohensD,Group_labels,'Color',cl_mat,'Sorting',false)
% boxplot_table_sorted_by_median(T,Group_labels,'SortMode','Group')

[figOut, ord] = plotCohensD_byGroup(cohensD, Group_labels, ...
    ...         % tells function to plot into this subplot
    'Colors', cl_mat, ...        % supply full color matrix matching canonical labels
    'Sorting', true, ...
    'ShowLegend', true,...
    'ShowXlabel', true ,...
    'GroupGap',2,...
    'FigureTitle',labels{1});
        
%% Clustering
% Apply parcellation
atlasname = 'canlab2024_coarse_fmriprep20_2mm';
atlas = load_atlas(atlasname);
f_correct = cell(1,30);
[contrastCells] = apply_parc_contrasiobj(data_cell,atlas,f_correct);  

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

[tbl, X, groupVec, varNames] = createMixedDesign(P_labels, cohensD_subjerr.n, cohensD_subjerr.d, {'mean_corr'});

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