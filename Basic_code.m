addpath(genpath('C:\Work\Toolboxes'))
addpath('C:\Work\NPS_metastudy\common_functions')
addpath('C:\Work\Toolboxes_general\spm12')
addpath('\\dartfs-hpc.dartmouth.edu\rc\lab\C\CANlab\labdata\projects\canlab_single_trials_for_git_repo')

%% Behavioral Datasets
% Non Placebo studies
ColNames_mat = load('C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\Postprocessing\labels_update8.mat');
ColNames = ColNames_mat.ColNames;
Group_labels = ColNames_mat.Group_labels;

group_ = cellfun(@(x) strcmp(x,'Placebo')||strcmp(x,'Placebo+C'), Group_labels);
group_id = double(group_)+1;
pl_labels = {'Not Placebo'; 'Placebo'};
pl_labels_list = pl_labels(group_id);

studydir = 'C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\included_studies';
studydir2 = 'C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\zunhammer2018_studies';

% Behavioral
% ccode = {[1],[3],[4],[21],[1],[11,12],[1],[1:2],[15:17],[9],[19],[1],[3],[1:2],[1],[1],[1],[1],[1]}';
ccode = {[1],[3],[4],[21],[1],[12],[1:2],[1],[1:2],[17],[9],[19],[1],[3],[1:2],[1],[1],[1],[1],[1],[1]}';

% Imaging
% idx_vec = [-1 1 1 1 -1 -1 1 1 1 1 -1 -1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1]'; % 1 = low pain - high pain
idx_vec = [-1 1 1 1 -1 1 1 1 1 1 1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1 1]';    % 1 = low pain - high pain

% Idx = {[1], [1], [1], [2], [1], [5:6], [1], [1:2], [1:3], [1],[1], [1], [1], [1:2], [1],[2],[1],[1],[1]}';
Idx = {[1], [1], [1], [2], [1], [6], [1:2], [1], [1:2], [3], [1],[1], [1], [1], [1:2], [1],[2],[1],[1],[1],[1]}';

%% Behavioral
[BigDat, ~] = NPSMS_harvest_canlab_cols(studydir, ccode, 55);

% Placebo studies
ccode2 = num2cell(ones(1,16));
[BigDat2, ~] = NPSMS_harvest_canlab_cols(studydir2, ccode2, 461);
BigDat = cat(2,BigDat, BigDat2);
lastRow = find(any(~isnan(BigDat),2), 1, 'last');
if ~isempty(lastRow); BigDat = BigDat(1:lastRow, :);end
% BigDat = l2normalize_columns(BigDat,'L2');
BigDat = l2normalize_columns(BigDat,'InvStd');

% Switch signs
BigDat_n = BigDat;
idx_vec2 = sign(mean(BigDat_n,'omitmissing'));
idx_mat = repmat(idx_vec2,[size(BigDat,1),1]);
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

% % Renorm
% std24 = BigDat_f{24};
% tmpl = image_vector('image_names', 'fmriprep20_template.nii');
% std24 = resample_space(tmpl,std24);
% BigDat_f{24} = std24; 

% Normalize
data_cell = alignFmriDataToReference(BigDat_f, [2]);

% % Num subjects
% data_cell2 = data_cell;
% data_cell2([5 8 11 18])=[];
% ns = cellfun(@(x) size(x.dat,2),data_cell2)

%% Preprocessing/Cleanup
% With sub-grid jitter to break the discrete rank grid
% V_out = harmonize_zero_preserve(V_in, 'jitter', 0.25);
% dat_obj = rescale(dat_obj, 'zeropreservequantile', 'jitter', 0.25, 'seed', 42);
% Harmonize and negate
% data_cell_rn = cellfun(@(x) rescale(x,'zeropreservequantile', 'jitter', 0.25,'seed', 42),data_cell,'UniformOutput',false);

fun = @(X) normalize_gm_by_wm_csf(X);
data_cell2 = cellfun(fun, data_cell, 'UniformOutput', false); % CSF removal on unnormalized data

% data_cell_rn2 = cellfun(@(x) harmonize_zero_preserve(x,'jitter', 0.25), data_cell, 'UniformOutput', false);  % CSF removal on normalized data
% fun_2 = @(X) normalize_gm_by_wm_csf(X,'do_scale',false);
% data_cell3 = cellfun(fun_2, data_cell, 'UniformOutput', false); % CSF removal on unnormalized data
% % fig = create_wm_fig(data_cell,data_cell2,data_cell3,Group_labels);
% [fig, stats] = create_wm_fig_cellfits(data_cell, data_cell2, data_cell3, Group_labels);

% GLM
[tmap_iv, pmap_iv, df,ts] = voxelwiseLM(data_cell2, C, st_vec, 1);
tmpl = image_vector('image_names', 'fmriprep20_template.nii');
tiv_rs = resample_space(tmap_iv, tmpl);

mapmin = min(tiv_rs.dat(:), [], 'omitnan');
mapmax = max(tiv_rs.dat(:), [], 'omitnan');
dat_disp = threshold(tiv_rs, [ts(2) ts(1)], 'raw-outside');
o2 = montage(dat_disp, 'trans', 'cmaprange', [mapmin ts(2) ts(1) mapmax]);
% montage(tiv_rs,'trans');

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

%% Studywise NPS-SIIPS analyses
ord = load('C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\Postprocessing\ratings_sort.mat');
ord = ord.ord;
Sgn_dat = apply_all_signatures(data_cell2,'conditionnames',ColNames);
% Sgn_dat2 = apply_nps(data_cell,'conditionnames',ColNames);

cl_mat = [0	113	188
216	82	24
236	176	31
125	46	141
118	171	47
76	189	237
161	19	46
247 129 191]/255;

cl_mat = cl_mat(1:length(unique(Group_labels)),:);

% labels = {'NPS','NPSpos','NPSneg','SIIPS','PINES','Rejection','VPS','FM_Multisens','FM_pain','Empathic_Care'};
labels = {'NPS','SIIPS'};
% ax = figure('Color','w','Position',[100 100 2400 1200]);
% hold on
% eval(sprintf('cohensD = cohensD_table_wrapper(Sgn_dat.%s,);',labels{1}))
cohensD = cohensD_table_wrapper(Sgn_dat.NPS,'Center','mean');

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
    'GroupRegions', true, 'Color',cl_mat,'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLim', [-2.25 1], 'YLim', [-2.25 1],'XLabel','NPS','YLabel','SIIPS');

% Circstat
% results_cov = classify_activation_pattern(-cohensD.d, -cohensD2.d, Group_labels);
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
    'FigureTitle','Pain score');

% NPS vs Pain
plotCohensD_bivariate_byGroup(cohensD,cohensD_rawpain, Group_labels, ...
    'GroupRegions', true,  'Color',cl_mat,'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLim', [-2.75 1], 'YLim', [-3 0.5],'XLabel','NPS','YLabel','Pain');

% SIIPS vs Pain
plotCohensD_bivariate_byGroup(cohensD2,cohensD_rawpain, Group_labels, ...
    'GroupRegions', true,  'Color',cl_mat,'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLim', [-2.75 1], 'YLim', [-3 0.5],'XLabel','SIIPS','YLabel','Pain');

results = classify_activation_pattern_v2( ...
    cohensD.d, cohensD.SE, ...
    cohensD2.d, cohensD2.SE, ...
    Group_labels,'color',cl_mat);
%% Subject-wise analyses
cl_mat2 = cat(1,[0.3 0.3 0.3],cl_mat);
cl_mat2(6,:) = [];

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

% [tbl, X, groupVec, varNames] = createMixedDesignEffectTable(Group_labels, cohensD.n,[],[],'referenceLevel','Placebo');
% [tbl, X, groupVec, varNames, info] = createMixedDesignDummyTable(Group_labels, cohensD.n,[],[],'referenceLevel','Placebo');

[tbl, X, groupVec, varNames, info] = createMixedDesignFullDummyTable(Group_labels, cohensD.n);

% Basic NPS
tbl.NPS = zscore(NPS_sc);
lme2 = fitlme(tbl, 'NPS ~ -1+ Conditioning + Mindfulness + Control + Placebo + Cognitive + Placebo_C + Remifentanil + Social + PCalib + (1|StudyID)');
disp(lme2);
[~,~,stats] = fixedEffects(lme2,'dfmethod','satterthwaite')
plotFixedEffects(stats(1:8,:), cl_mat,'specialIntercept', false)

% 
% X = designMatrix(lme2, 'Fixed');


% Basic SIIPS
tbl.SIIPS = zscore(SIIPS_sc);
lme2 = fitlme(tbl, 'SIIPS  ~ -1+ Conditioning + Mindfulness + Control + Cognitive + Placebo_C + Remifentanil + Social + Placebo  + PSite + (1|StudyID)');
disp(lme2);
[~,~,stats] = fixedEffects(lme2,'dfmethod','satterthwaite')
plotFixedEffects(stats(1:8,:), cl_mat,'specialIntercept', false)

% NPS~SIIPS 
lme2 = fitlme(tbl, ['NPS ~ -1+ Conditioning + Mindfulness + Control + Cognitive + Placebo_C + Remifentanil + Social + Placebo +' ...
    'SIIPS:Conditioning + SIIPS:Mindfulness + SIIPS:Control + SIIPS:Cognitive + SIIPS:Placebo_C + SIIPS:Remifentanil + SIIPS:Social + SIIPS:Placebo + (1|StudyID)']);
disp(lme2);
[~,~,stats] = fixedEffects(lme2,'dfmethod','satterthwaite')
plotFixedEffects(stats(9:end,:), cl_mat,'specialIntercept', false)

% % Regressed regression with JKCorr
% tbl.NPS = (NPS_sc_reg);
% lme2 = fitlme(tbl, 'NPS ~ 1+ Label_Conditioning + Label_Mindfulness + Label_Cognitive + Label_Placebo_ + Label_Remifentanil + Label_Social + (1|StudyID)');
% % lme2 = fitlme(tbl, 'NPS ~ 1+ Label_Conditioning + Label_Mindfulness + Label_Cognitive + Label_Placebo_ + Label_Remifentanil + Label_Social + (1|StudyID)','Weights',(sim_vect+1)/2);
% disp(lme2);
% [~,~,stats] = fixedEffects(lme2,'dfmethod','satterthwaite')

%% Subject-wise analyses - raw pain/NPS
% With raw pain
rem_st = 3;

pain_sc = BigDat_n;
pain_sc(:,rem_st) = [];
pain_sc = pain_sc(~isnan(pain_sc));

% Subjectwise NPS
NPS_scr = Sgn_dat.NPS;
NPS_scr(:,rem_st)=[];
NPS_scr = NPS_scr{:,:};
NPS_scr = NPS_scr(:);
NPS_scr(isnan(NPS_scr))=[];

% Subjectwise SIIPS
SIIPS_scr = Sgn_dat.SIIPS;
SIIPS_scr(:,rem_st)=[];
SIIPS_scr = SIIPS_scr{:,:};
SIIPS_scr = SIIPS_scr(:);
SIIPS_scr(isnan(SIIPS_scr))=[];

% [tbl, X, groupVec, varNames] = createMixedDesignEffectTable(Group_labels, cohensD.n,[],[],'referenceLevel','Placebo');
grp2 = Group_labels;
grp2(rem_st)=[];
num2 = cohensD.n;
num2(rem_st)=[];
% [tbl2,~, groupVec2, ~, ~] = createMixedDesignDummyTable(grp2, num2,[],[],'referenceLevel','Placebo');
[tbl2, X, groupVec, varNames, info] = createMixedDesignFullDummyTable(grp2, num2);

% Fixed effect
tbl2.pain = pain_sc;
tbl2.NPS = zscore(NPS_scr);
tbl2.SIIPS = zscore(SIIPS_scr)
lme2 = fitlme(tbl2, 'pain ~ -1+ NPS + SIIPS + (1|StudyID)');
disp(lme2);
[~,~,stats] = fixedEffects(lme2,'dfmethod','satterthwaite')
cl = lines(3);
plotFixedEffects(stats, cl,'specialIntercept', false)

lme2 = fitlme(tbl2, 'pain ~ -1+ Conditioning + Mindfulness + Control + Cognitive + Placebo_C + Remifentanil + Social + Placebo + PCalib+ (1|StudyID)');
disp(lme2);
[~,~,stats] = fixedEffects(lme2,'dfmethod','satterthwaite')
plotFixedEffects(stats(1:8,:), cl_mat,'specialIntercept', false)

lme2 = fitlme(tbl2, ['pain ~ -1+ Conditioning + Mindfulness + Control + Cognitive + Placebo_C + Remifentanil + Social + Placebo +' ...
    'NPS:Conditioning + NPS:Mindfulness + NPS:Control + NPS:Cognitive + NPS:Placebo_C + NPS:Remifentanil + NPS:Social + NPS:Placebo +'...
    'SIIPS:Conditioning + SIIPS:Mindfulness + SIIPS:Control + SIIPS:Cognitive + SIIPS:Placebo_C + SIIPS:Remifentanil + SIIPS:Social + SIIPS:Placebo + (1|StudyID)']);
disp(lme2);
[~,~,stats] = fixedEffects(lme2,'dfmethod','satterthwaite')
% NPS
plotFixedEffects(stats(9:16,:), cl_mat,'specialIntercept', false)
% SIIPS
plotFixedEffects(stats(17:end,:), cl_mat,'specialIntercept', false)

% % PosvNeg
% plotCohensD_bivariate_byGroup(cohensD_pos,cohensD_neg, Group_labels, ...
%     'GroupRegions', true, 'color',cl_mat,'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLabel','Pos Axis','YLabel','Neg Axis');

CP = tbl2.NPS + 1i*tbl2.SIIPS;
tbl2.CP_mag = abs(CP);
tbl2.CP_tan = angle(CP);

lme2 = fitlme(tbl2, 'CP_tan ~ -1+ Conditioning + Mindfulness + Control + Cognitive + Placebo_C + Remifentanil + Social + Placebo + (1|StudyID)');
disp(lme2);
[~,~,stats] = fixedEffects(lme2,'dfmethod','satterthwaite')
plotFixedEffects(stats, cl_mat,'specialIntercept', false)


%% Common factors +ve -ve

% Apply parcellation
atlasname = 'canlab2024_coarse_fmriprep20_2mm';
atlas = load_atlas(atlasname);
f_correct = cell(1,length(data_cell2));
[contrastCells] = apply_parc_contrasiobj(data_cell2,atlas,f_correct);  
nstud = length(data_cell2);
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
 
% % C 3 is pos
% idx_pos2 = false(1,length(idx_n));
% idx_pos2(idx(stats.all_cluster_labels(:,2)==3)) = true;
% idx_neg2 = false(1,length(idx_n));
% idx_neg2(idx(or(stats.all_cluster_labels(:,2)==1,stats.all_cluster_labels(:,2)==2))) = true;
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
    'GroupRegions', true, 'color',cl_mat,'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLabel','Pos Axis','YLabel','Neg Axis');

% % cohensD_neg.d = -(median(studymean_reg_(:,idx_neg),2));
% plotCohensD_bivariate_byGroup(cohensD_pos,cohensD_neg, Group_labels, ...
%     'GroupRegions', true, 'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLabel','Pos','YLabel','Neg');

% Circstats
% results = classify_activation_pattern(cohensD_pos.d, cohensD_neg.d, Group_labels);
results = classify_activation_pattern_v2( ...
    cohensD_pos.d, cohensD_pos.SE, ...
    cohensD_neg.d, cohensD_neg.SE, ...
    Group_labels,'color',cl_mat);
[d,p,stats2] = manova1([cohensD_pos.d, cohensD_neg.d],Group_labels);

% postStats = manovaPairwisePostHoc(Y, Group_labels, 'Adjust', 'holm');
% cl_mat = lines(numel(postStats.Estimate));
% plotFixedEffects(postStats, cl_mat);

% Circstats - placebo
results = classify_activation_pattern_v2( ...
    cohensD_pos.d, cohensD_pos.SE, ...
    cohensD_neg.d, cohensD_neg.SE, ...
    pl_labels_list,'color',cl_mat);

% Corr with raw pain?
plotCohensD_bivariate_byGroup(cohensD_pos,cohensD_rawpain, Group_labels, ...
    'GroupRegions', true, 'color',cl_mat,'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLim', [-0.3 0.6], 'YLim', [-3 0.2], 'XLabel','Pos Axis','YLabel','Pain Score');

% Corr with raw pain?
plotCohensD_bivariate_byGroup(cohensD_neg,cohensD_rawpain, Group_labels, ...
    'GroupRegions', true, 'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLim', [-0.3 0.6], 'YLim', [-3 0], 'XLabel','Neg Axis','YLabel','Pain Score');


% % Cluster map
% idx_vec = zeros(1,length(idx_n));
% kmax = 3;
% for ii = 1:kmax
% idx_vec(idx(stats.all_cluster_labels(:,kmax-1)==ii)) =ii;
% end
% 
% idx_vec(idx_vec==3)=0;
% idx_vec(idx_vec==1)=3;
% 
% atlas    = load_atlas('canlab2024_coarse_fmriprep20_2mm');
% tiv_rs   = roi_vector_to_image_vector(atlas, idx_vec);
% plot_signed_threshold(tiv_rs, 0.1);
% 
% montage(tiv_rs)

%% Subject Wise
setid = setdiff(1:length(contrastCells),rem_st);
CMain = vertcat(contrastCells{setid});
CPos = CMain(:,idx_pos2);
CNeg = CMain(:,idx_neg2);

[tbl2,~, groupVec2, ~, ~] = createMixedDesignFullDummyTable(grp2, num2);
% Fixed effect
tbl2.pain = zscore(pain_sc);
tbl2.NPS = zscore(NPS_scr);
tbl2.SIIPS = zscore(SIIPS_scr);
tbl2.CPos = zscore(mean(CPos,2));
tbl2.CNeg = zscore(-mean(CNeg,2));

lme_grp = fitlme(tbl2, ...
    ['CNeg ~ -1 +' ...
        'Conditioning + Mindfulness + Control + ' ...
        'Cognitive + Placebo_C + Remifentanil + Social + Placebo +' ...
     '(1|StudyID)']);
disp(lme_grp);
[~,~,stats] = fixedEffects(lme_grp,'dfmethod','satterthwaite');
plotFixedEffects(stats,  cl_mat,'specialIntercept', false)


lme_grp = fitlme(tbl2, ...
    ['CPos ~ -1 +' ...
        'Conditioning + Mindfulness + Control + ' ...
        'Cognitive + Placebo_C + Remifentanil + Social + Placebo +' ...
     '(1|StudyID)']);

disp(lme_grp);
[~,~,stats] = fixedEffects(lme_grp,'dfmethod','satterthwaite');
plotFixedEffects(stats,  cl_mat,'specialIntercept', false)

% lme_grp = fitlme(tbl2, ...
%     ['CPos ~ -1 +' ...
%          'CNeg*Conditioning + CNeg*Mindfulness + CNeg*Control + ' ...
%         'CNeg*Cognitive + CNeg*Placebo_C + CNeg*Remifentanil + CNeg*Social + CNeg*Placebo +' ...

lme_grp = fitlme(tbl2, ['CPos ~ -1+ Conditioning + Mindfulness + Control + Cognitive + Placebo_C + Remifentanil + Social + Placebo +' ...
    'CNeg:Conditioning + CNeg:Mindfulness + CNeg:Control + CNeg:Cognitive + CNeg:Placebo_C + CNeg:Remifentanil + CNeg:Social + CNeg:Placebo + (1|StudyID)']);
disp(lme_grp);
[~,~,stats] = fixedEffects(lme_grp,'dfmethod','satterthwaite');
plotFixedEffects(stats(9:end,:),  cl_mat,'specialIntercept', false)

lme_grp = fitlme(tbl2, ...
    ['CPos ~ 1 + CNeg +' ...
     '(1|StudyID)']);
disp(lme_grp);
[~,~,stats] = fixedEffects(lme_grp,'dfmethod','satterthwaite');
plotFixedEffects(stats(9:end,:), cl_mat2 )

lme_grp = fitlme(tbl2, ...
    ['pain ~ -1 + ' ...
     'Conditioning + Mindfulness + Control + ' ...
        'Cognitive + Placebo_C + Remifentanil + Social + ' ...
    'CPos:Conditioning + CPos:Mindfulness + CPos:Control + ' ...
        'CPos:Cognitive + CPos:Placebo_C + CPos:Remifentanil + CPos:Social + ' ...
    'CNeg:Conditioning + CNeg:Mindfulness + CNeg:Control + ' ...
        'CNeg:Cognitive + CNeg:Placebo_C + CNeg:Remifentanil + CNeg:Social + ' ...
     '(1|StudyID)']);

[~,~,stats] = fixedEffects(lme_grp,'dfmethod','satterthwaite');
% POS
plotFixedEffects(stats([8:14],:),  cl_mat,'specialIntercept', false)
% NEG
plotFixedEffects(stats([15:end],:),  cl_mat,'specialIntercept', false)

lme2 = fitlme(tbl2, 'pain ~ 1+ CPos + CNeg + (1|StudyID)');
disp(lme2);
[~,~,stats] = fixedEffects(lme2,'dfmethod','satterthwaite')
cl = lines(3);
plotFixedEffects(stats, cl)

lme2 = fitlme(tbl2, 'SIIPS ~ 1+ CPos + CNeg + (1|StudyID)');
disp(lme2);
[~,~,stats] = fixedEffects(lme2,'dfmethod','satterthwaite')
cl = lines(3);
plotFixedEffects(stats, cl)

CP = tbl2.CPos + 1i*tbl2.CNeg;
CP = tbl2.NPS + 1i*tbl2.SIIPS;

tbl2.CP_mag = abs(CP);
tbl2.CP_tan = angle(CP);

lme_grp = fitlme(tbl2, ...
    ['CP_mag ~ -1 +' ...
        'Conditioning + Mindfulness + Control + ' ...
        'Cognitive + Placebo_C + Remifentanil + Social + Placebo +' ...
     '(1|StudyID)']);
disp(lme_grp);
[~,~,stats] = fixedEffects(lme_grp,'dfmethod','satterthwaite');
plotFixedEffects(stats, cl_mat,'specialIntercept', false)

%% Neural ROIs
% Apply parcellation
atlasname = 'canlab2024_coarse_fmriprep20_2mm';
atlas = load_atlas(atlasname);
f_correct = cell(1,length(data_cell2));
[contrastCells] = apply_parc_contrasiobj(data_cell2,atlas,f_correct);  
nstud = length(data_cell2);
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

% idx_1 = idx_pos;
% idx_2 = idx_neg;
 
% C 3 is pos
idx_1 = false(1,length(idx_n));
idx_1(idx(stats.all_cluster_labels(:,2)==1)) = true;
idx_2 = false(1,length(idx_n));
idx_2(idx(stats.all_cluster_labels(:,2)==2)) = true;
idx_3 = false(1,length(idx_n));
idx_3(idx(stats.all_cluster_labels(:,2)==3)) = true; % Positive

% llist_pos = atlas.label_descriptions(idx_pos);
% sum(ismember(atlas.label_descriptions(idx_1),llist_pos))

idx_pos2 = idx_3;
% idx_neg2 = or(idx_1,idx_2);

idx_pos2 = idx_1;
idx_neg2 = or(idx_2,idx_2);

cohensD_pos = cohensD;
cohensD_pos.d = -median(studymean_reg_(:,idx_pos2),2);
cohensD_pos.SE = median(studysem_reg_(:,idx_pos2),2);

cohensD_neg = cohensD;
cohensD_neg.d = -median(studymean_reg_(:,idx_neg2),2);
cohensD_neg.SE = median(studysem_reg_(:,idx_neg2),2);

% PosvNeg
plotCohensD_bivariate_byGroup(cohensD_pos,cohensD_neg, Group_labels, ...
    'GroupRegions', true, 'color',cl_mat,'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLabel','Pos Axis','YLabel','Neg Axis');

% % cohensD_neg.d = -(median(studymean_reg_(:,idx_neg),2));
% plotCohensD_bivariate_byGroup(cohensD_pos,cohensD_neg, Group_labels, ...
%     'GroupRegions', true, 'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLabel','Pos','YLabel','Neg');

% Circstats
% results = classify_activation_pattern(cohensD_pos.d, cohensD_neg.d, Group_labels);
results = classify_activation_pattern_v2( ...
    cohensD_pos.d, cohensD_pos.SE, ...
    cohensD_neg.d, cohensD_neg.SE, ...
    Group_labels,'color',cl_mat);
[d,p,stats2] = manova1([cohensD_pos.d, cohensD_neg.d],Group_labels);

% postStats = manovaPairwisePostHoc(Y, Group_labels, 'Adjust', 'holm');
% cl_mat = lines(numel(postStats.Estimate));
% plotFixedEffects(postStats, cl_mat);

% Circstats - placebo
results = classify_activation_pattern_v2( ...
    cohensD_pos.d, cohensD_pos.SE, ...
    cohensD_neg.d, cohensD_neg.SE, ...
    pl_labels_list,'color',cl_mat);

% Corr with raw pain?
plotCohensD_bivariate_byGroup(cohensD_pos,cohensD_rawpain, Group_labels, ...
    'GroupRegions', true, 'color',cl_mat,'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLim', [-0.3 0.6], 'YLim', [-3 0.2], 'XLabel','Pos Axis','YLabel','Pain Score');

% Corr with raw pain?
plotCohensD_bivariate_byGroup(cohensD_neg,cohensD_rawpain, Group_labels, ...
    'GroupRegions', true, 'GroupRegionStyle', 'ellipse','GroupRegionAlpha', 0.08,'XLim', [-0.3 0.6], 'YLim', [-3 0], 'XLabel','Neg Axis','YLabel','Pain Score');

setid = setdiff(1:length(contrastCells),rem_st);
CMain = vertcat(contrastCells{setid});
C1 = CMain(:,idx_1);
C2 = CMain(:,idx_2);
C3 = CMain(:,idx_3);

tbl2.C1 = zscore(mean(C1,2));
tbl2.C2 = zscore(mean(C2,2));
tbl2.C3 = zscore(mean(C3,2));
% % Cluster map
% idx_vec = zeros(1,length(idx_n));
% kmax = 3;
% for ii = 1:kmax
% idx_vec(idx(stats.all_cluster_labels(:,kmax-1)==ii)) =ii;
% end
% 
% idx_vec(idx_vec==3)=0;
% idx_vec(idx_vec==1)=3;
% 
% atlas    = load_atlas('canlab2024_coarse_fmriprep20_2mm');
% tiv_rs   = roi_vector_to_image_vector(atlas, idx_vec);
% plot_signed_threshold(tiv_rs, 0.1);
% 
% montage(tiv_rs)

%% Decoding analyses
% Bar 2: same, averaged within each intervention (all pairs containing it).
% Dashed red line = isotropy reference b_iso = (k-1)*a.
addpath('common_functions');
S = load('datamat.mat'); cl_mat = S.cl_mat;

R = run_test1('MODE','all','COORDS',{'C1','C2','C3'},'N_PERM',2000,'N_BOOT',1000);

plotTest1AB(R, cl_mat, 'Test 1 - on-axis vs null-space separation (all pairs)');
exportgraphics(gcf, 'fig_test1_ab_overall.png', 'Resolution', 200);

plotTest1AB(R, cl_mat, 'Test 1 - separation split by intervention', 'byGroup', true);
exportgraphics(gcf, 'fig_test1_ab_bygroup.png', 'Resolution', 200);

disp('test1 figures written');
