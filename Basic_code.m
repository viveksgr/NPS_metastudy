addpath(genpath('C:\Work\Toolboxes'))
addpath('C:\Work\NPS_metastudy\common_functions')
addpath('C:\Work\Toolboxes_general\spm12')
addpath('\\dartfs-hpc.dartmouth.edu\rc\lab\C\CANlab\labdata\projects\canlab_single_trials_for_git_repo')

%% Behavioral Datasets
% Non Placebo studies
ColNames_mat = load('C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\Postprocessing\labels_update3.mat');
ColNames = ColNames_mat.ColNames;
Group_labels = ColNames_mat.Group_labels;

group_ = cellfun(@(x) strcmp(x,'Placebo')||strcmp(x,'Placebo+'), Group_labels);
group_id = double(group_)+1;
pl_labels = {'Not Placebo'; 'Placebo'};
pl_labels_list = pl_labels(group_id);

studydir = 'C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\included_studies';
studydir2 = 'C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\zunhammer2018_studies';

%% Behavioral
ccode = {[1],[3],[4],[21],[1],[11,12],[15:17],[9],[19],[1],[3],[1:2],[1]}';
[BigDat, ~] = NPSMS_harvest_canlab_cols(studydir, ccode, 55);

% Placebo studies
ccode2 = num2cell(ones(1,16));
[BigDat2, ~] = NPSMS_harvest_canlab_cols(studydir2, ccode2, 55);
BigDat = cat(2,BigDat, BigDat2);
lastRow = find(any(~isnan(BigDat),2), 1, 'last');
if ~isempty(lastRow); BigDat = BigDat(1:lastRow, :);end
BigDat = l2normalize_columns(BigDat,'L2');
% BigDat = l2normalize_columns(BigDat,'InvStd');

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

% For ILCP
% idx_vec = [-1 1 1 1 -1 -1 1 1 1 -1 -1 -1 -1 1 1 1 -1 -1 -1]'; % 1 = low pain - high pain
% Idx = {[1], [1], [1], [2], [1], [5:6], [1:2], [1:3], [1],[1], [1], [1], [1:2], [2]}';

idx_vec = [-1 1 1 1 -1 -1 1 -1 -1 -1 -1 1 1 1 -1 -1 -1]'; % 1 = low pain - high pain
Idx = {[1], [1], [1], [2], [1], [5:6], [1:3], [1],[1], [1], [1], [1:2], [2]}';

st_vec_ind = cellfun(@(x) numel(x),Idx);
st_vec = [];
kk = 0;
for ii=1:length(st_vec_ind) 
    for jj = 1:st_vec_ind(ii)
        kk = kk+1;
        st_vec(kk)=ii;
    end
end
st_vec = cat(2,st_vec,st_vec(end)+1:st_vec(end)+16);

% Dummy contrast
C = group_id;

[BigDat_f] = harvest_canlab_funccols(studydir, Idx);
ncontrast = numel(BigDat_f);
for ii = 1:ncontrast
    % BigDat_f{ii}=  rescale(BigDat_f{ii}, 'csf_mean_var');
    BigDat_f{ii}=  rescale(BigDat_f{ii}, 'l2norm_images');
    % BigDat_f{ii}=  rescale(BigDat_f{ii}, 'prctileimages');
    BigDat_f{ii}.dat  = BigDat_f{ii}.dat.*idx_vec(ii); % Change signs of the contrasts so that all contrasts are for analgesia (low pain - high pain)
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
% With sub-grid jitter to break the discrete rank grid
% V_out = harmonize_zero_preserve(V_in, 'jitter', 0.25);

% dat_obj = rescale(dat_obj, 'zeropreservequantile', 'jitter', 0.25, 'seed', 42);

% Harmonize and negate
% data_cell_rn = cellfun(@(x) rescale(x,'zeropreservequantile', 'jitter', 0.25,'seed', 42),data_cell,'UniformOutput',false);

fun = @(X) normalize_gm_by_wm_csf(X);
data_cell2 = cellfun(fun, data_cell, 'UniformOutput', false); % CSF removal on unnormalized data
data_cell_rn2 = cellfun(@(x) harmonize_zero_preserve(x,'jitter', 0.25), data_cell, 'UniformOutput', false);  % CSF removal on normalized data
fun_2 = @(X) normalize_gm_by_wm_csf(X,'do_scale',false);
data_cell3 = cellfun(fun_2, data_cell, 'UniformOutput', false); % CSF removal on unnormalized data

% fig = create_wm_fig(data_cell,data_cell2,data_cell3,Group_labels);
[fig, stats] = create_wm_fig_cellfits(data_cell, data_cell2, data_cell3, Group_labels);


% GLM
[tmap_iv, pmap_iv, df,ts] = voxelwiseLM(data_cell_rn2, C, st_vec, 2);
% montage(tmap_iv)
% GLM
% [tmap_iv2, pmap_iv2, df,ts] = voxelwiseLM(data_cell_rn, C, st_vec, 1);
% GLM-2
% [fig, stats] = scatter_density_unity(tmap_iv2.dat,tmap_iv.dat, 'xlabel', 'Before JK Corr', 'ylabel', 'After JK Corr');

% 
% % Normalization
% k = 0;
% figure('Position',[0.5 0.5 1280 720])
% hold on
% for ii = 1:30
%     k = k+1;
%     subplot(5,6,k)
%     plotGlobalMeans(data_cell_rn{ii})
%     ylim([-0.5 0.5])
%     title(ColNames{ii},'Interpreter','none')
% end

%% NPS masks
ord = load('C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\Postprocessing\ratings_sort.mat');
ord = ord.ord;
Sgn_dat = apply_all_signatures(data_cell2,'conditionnames',ColNames);
% Sgn_dat2 = apply_nps(data_cell,'conditionnames',ColNames);

cl_mat = lines(length(unique(Group_labels)));
% labels = {'NPS','NPSpos','NPSneg','SIIPS','PINES','Rejection','VPS','FM_Multisens','FM_pain','Empathic_Care'};
labels = {'NPS','SIIPS'};
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
    'GroupGap',1,...
    'FigureTitle',labels{1});

eval(sprintf('cohensD2 = cohensD_table_wrapper(Sgn_dat.%s);',labels{2}))
[figOut, ~] = plotCohensD_byGroup(cohensD2(ord,:), Group_labels(ord), ...
    ...         % tells function to plot into this subplot
    'Colors', cl_mat, ...        % supply full color matrix matching canonical labels
    'Sorting', true, ...
    'ShowLegend', true,...
    'ShowXlabel', true ,...
    'GroupGap',1,...
    'FigureTitle',labels{2});

% NPS vs SIIPS
plotCohensD_bivariate_byGroup(cohensD, cohensD2,Group_labels, ...
    'GroupRegions', true, 'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLim', [-1.5 1], 'YLim', [-1.5 1],'XLabel','NPS','YLabel','SIIPS');

% Circstat
results_cov = classify_activation_pattern(-cohensD.d, -cohensD2.d, Group_labels);
% [d,p,stats] = manova1([-cohensD.d, -cohensD2.d],Group_labels);


% Pain ratings
cohensD_rawpain = cohensD_table_wrapper(BigDat_n);
cohensD_rawpain.Name = cohensD.Name;
[figOut, ord] = plotCohensD_byGroup(cohensD_rawpain, Group_labels, ...
    ...         % tells function to plot into this subplot
    'Colors', cl_mat, ...        % supply full color matrix matching canonical labels
    'Sorting', true, ...
    'ShowLegend', true,...
    'ShowXlabel', true ,...
    'GroupGap',1,...
    'FigureTitle',labels{1});


% NPS vs Pain
plotCohensD_bivariate_byGroup(cohensD,cohensD_rawpain, Group_labels, ...
    'GroupRegions', true, 'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLim', [-2.75 1], 'YLim', [-2.75 1],'XLabel','NPS','YLabel','Pain');

% SIIPS vs Pain
plotCohensD_bivariate_byGroup(cohensD2,cohensD_rawpain, Group_labels, ...
    'GroupRegions', true, 'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLim', [-2.75 1], 'YLim', [-2.75 1],'XLabel','SIIPS','YLabel','Pain');

% Subjectwise NPS
NPS_sc = Sgn_dat.NPS;
NPS_sc = NPS_sc{:,:};
NPS_cell = conv2cell(NPS_sc); 
NPS_sc = NPS_sc(:);
NPS_sc(isnan(NPS_sc))=[];

SIIPS_sc = Sgn_dat.SIIPS;
SIIPS_sc = SIIPS_sc{:,:};
SIIPS_sc = SIIPS_sc(:);
SIIPS_sc(isnan(SIIPS_sc))=[];
%% Jackknife similarity
% [sim_values, d, low_agreement, Nvox] = jackknife_similarity(data_cell_rn{1});
[sim_values, ~, d] = cellfun(@(x) jackknife_similarity(x),data_cell_rn2,'UniformOutput',false);
[sim_values_sad, ~, d1] = cellfun(@(x) jackknife_similarity(x,'similarity_metric','standardized_abs_deviation'),data_cell_rn2,'UniformOutput',false);
[sim_values_msz, ~, d2] = cellfun(@(x) jackknife_similarity(x,'similarity_metric','mean_shift_z'),data_cell_rn2,'UniformOutput',false);
[sim_values_ssz, ~, d3] = cellfun(@(x) jackknife_similarity(x,'similarity_metric','scale_shift_z'),data_cell_rn2,'UniformOutput',false);
% 
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

% [tbl, X, groupVec, varNames] = createMixedDesignEffectTable(Group_labels, cohe    nsD.n,[],[],'referenceLevel','Cognitive');
% [tbl, X, groupVec, varNames] = createMixedDesignEffectTable(Group_labels, cohensD.n,[],[],'referenceLevel','Placebo');
[tbl, X, groupVec, varNames, info] = createMixedDesignDummyTable(Group_labels, cohensD.n,[],[],'referenceLevel','Placebo');

% Basic NPS
tbl.NPS = (NPS_sc);
lme2 = fitlme(tbl, 'NPS ~ 1+ Label_Conditioning + Label_Mindfulness + Label_Cognitive + Label_Placebo_ + Label_Remifentanil + Label_Social + (1|StudyID)');
disp(lme2);
[~,~,stats] = fixedEffects(lme2,'dfmethod','satterthwaite')

% Regressed regression
tbl.NPS = (NPS_sc_reg);
% lme2 = fitlme(tbl, 'Y ~ 1+ Label_Conditioning + Label_Mindfulness + Label_Placebo + Label_Placebo_ + Label_Remifentanil + Label_Social + (1|StudyID)');
lme2 = fitlme(tbl, 'NPS ~ 1+ Label_Conditioning + Label_Mindfulness + Label_Cognitive + Label_Placebo_ + Label_Remifentanil + Label_Social + (1|StudyID)','Weights',(sim_vect+1)/2);
disp(lme2);
[~,~,stats] = fixedEffects(lme2,'dfmethod','satterthwaite')
cl_mat2 = cat(1,[0.3 0.3 0.3],cl_mat);
cl_mat2(5,:) = [];
plotFixedEffects(stats, cl_mat2)

% For rawpain
[tbl, X, groupVec, varNames] = createMixedDesignEffectTable(Group_labels, cohensD_rawpain.n,[],[],'referenceLevel','Placebo');

% NPS vs SIIPS
NPS_sc_reg = regress_out(NPS_sc, sim_vect);
SIIPS_sc_reg = regress_out(SIIPS_sc, sim_vect);

results_cov = classify_activation_pattern(-NPS_sc_reg, -SIIPS_sc_reg, groupVec);
%% Common factors +ve -ve

% Apply parcellation
atlasname = 'canlab2024_coarse_fmriprep20_2mm';
atlas = load_atlas(atlasname);
f_correct = cell(1,length(data_cell2));
[contrastCells] = apply_parc_contrasiobj(data_cell2,atlas,f_correct);  
nstud = length(data_cell_rn2);
c_map = extract_roi_averages(tmap_iv, atlas);
roi_means = zeros(1,length(atlas.labels));
for nn=1:length(atlas.labels); roi_means(nn)=median(c_map(nn).dat); end


idx_pos = roi_means>ts(1);
idx_neg = roi_means<-ts(2);
% idx_pos = roi_means>3.1;
% idx_neg = roi_means<-3.1;
idx_n = 1:length(atlas.labels);
idx = [idx_n(idx_pos) idx_n(idx_neg)];

studymeans = cellfun(@(x) mean(x),contrastCells,'UniformOutput',false);
studysem = cellfun(@(x) std(x)./sqrt(size(x,1)),contrastCells,'UniformOutput',false);

studymean_reg_ = vertcat(studymeans{:});
studysem_reg_ = vertcat(studysem{:});
studymean_reg = studymean_reg_(:,idx);

studymean_corr = corrcoef(studymean_reg);

% studywise perm
stats = clusterdata_permtest(studymean_corr, ...
    'k', 2:7, ...
    'distancemetric', 'correlation', ...
    'linkagemethod', 'average', ...
    'reducedims', true, ...
    'nperm', 100, ...
    'doplot', true, ...
    'verbose', true);

[lid,argsort] = sort(stats.best_cluster_labels);
idxsort = idx(argsort);
nanat = length(stats.best_cluster_labels);
figure('Position',[0.5 0.5 640 480])
imagesc(1:nanat,1:nanat,studymean_corr(argsort,argsort));
hold on
xticks(1:nanat)
yticks(1:nanat)
xticklabels(atlas.labels(idxsort))
yticklabels(atlas.labels(idxsort))

k_picker = 2;
% Crosstab analysis
idx_id = [ones(sum(idx_pos),1); 2*ones(sum(idx_neg),1)];
crosstab_cluster_overlap(idx_id,stats.all_cluster_labels(:,k_picker-1),  'name1', 'GLM contrast', 'name2', 'K-means')

idx_pos2 = idx_pos;
idx_neg2 = idx_neg;
% 
% idx_pos2 = false(1,length(idx_n));
% idx_pos2(idx(stats.all_cluster_labels(:,2)==2)) = true;
% idx_neg2 = false(1,length(idx_n));
% idx_neg2(idx(or(stats.all_cluster_labels(:,2)==1,stats.all_cluster_labels(:,2)==1))) = true;
% llist_1 = atlas.label_descriptions(idx_pos2);
% llist_2 = atlas.label_descriptions(idx_neg2);


cohensD_pos = cohensD;
cohensD_pos.d = median(studymean_reg_(:,idx_pos2),2);
cohensD_pos.SE = median(studysem_reg_(:,idx_pos2),2);

cohensD_neg = cohensD;
cohensD_neg.d = -median(studymean_reg_(:,idx_neg2),2);
cohensD_neg.SE = median(studysem_reg_(:,idx_neg2),2);

% PosvNeg
plotCohensD_bivariate_byGroup(cohensD_pos,cohensD_neg, Group_labels, ...
    'GroupRegions', true, 'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLabel','Neg Axis - 1','YLabel','Neg Axis - 2');

% % cohensD_neg.d = -(median(studymean_reg_(:,idx_neg),2));
% plotCohensD_bivariate_byGroup(cohensD_pos,cohensD_neg, Group_labels, ...
%     'GroupRegions', true, 'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLabel','Pos','YLabel','Neg');

% Circstats
results = classify_activation_pattern(cohensD_pos.d, cohensD_neg.d, Group_labels);
[d,p,stats] = manova1([cohensD_pos.d, cohensD_neg.d],Group_labels);

% postStats = manovaPairwisePostHoc(Y, Group_labels, 'Adjust', 'holm');
% cl_mat = lines(numel(postStats.Estimate));
% plotFixedEffects(postStats, cl_mat);

% Circstats - placebo
results = classify_activation_pattern(cohensD_pos.d, cohensD_neg.d, pl_labels_list);
% [d,p,stats] = manova1([results.theta_deg results.theta_rad],double(group_));

% Corr with raw pain?
plotCohensD_bivariate_byGroup(cohensD_neg,cohensD_rawpain, Group_labels, ...
    'GroupRegions', true, 'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLim', [-0.75 0.3], 'YLim', [-3 0], 'XLabel','Neg Axis','YLabel','Pain Score');

results = classify_activation_pattern(cohensD_pos.d, cohensD_rawpain.d, pl_labels_list);


% Cluster map
idx_vec = zeros(1,length(idx_n));
kmax = 3;
for ii = 1:kmax
idx_vec(idx(stats.all_cluster_labels(:,kmax-1)==ii)) =ii;
end

idx_vec(idx_vec==3)=0;
idx_vec(idx_vec==1)=3;

atlas    = load_atlas('canlab2024_coarse_fmriprep20_2mm');
tiv_rs   = roi_vector_to_image_vector(atlas, idx_vec);
plot_signed_threshold(tiv_rs, 0.1);

montage(tiv_rs)



















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