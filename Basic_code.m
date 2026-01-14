addpath(genpath('C:\Work\Toolboxes'))
addpath(genpath('C:\Work\mind\toolboxes'))
addpath('C:\Work\NPS_metastudy\common_functions')

%% Behavioral Datasets
% Non Placebo studies
ColNames_mat = load('C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\Postprocessing\labels.mat');
ColNames = ColNames_mat.ColNames;
Group_labels = ColNames_mat.group_var;

ccode = {[1],[3],[21],[10,11,12],[15:18],[9],[19],[1],[3],[1]}';
studydir = 'C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\included_studies';
[BigDat, ~] = NPSMS_harvest_canlab_cols(studydir, ccode, 50);

% Placebo studies
studydir2 = 'C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\zunhammer2018_studies';
ccode2 = num2cell(ones(1,15));
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
BigDat_n = BigDat.*idx_mat;
plot_sorted_bar_with_errors(BigDat_n, ColNames)

%% Imaging Data
%  desc = descriptives(DATA_OBJ_CON{1, 1}  , ['noverbose', 'plotcoverage']);
%  stats = clusterdata_permtest(DATA_OBJ_CON{1, 1}.dat, 'k', [2:20], 'reducedims', true, 'ndims', 25);
%  % desc = descriptives(DATA_OBJ_CON{1, 1}  , ['noverbose', 'plotcoverage']);
%  % qc_metrics_second_level(DATA_OBJ_CON{1, 1});
% % help robfit_parcelwise
ColNames_mat = load('C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\Postprocessing\labels.mat');
ColNames = ColNames_mat.ColNames;
Group_labels = ColNames_mat.group_var;
idx_vec = [1 -1 -1 -1 1 -1 1 1 1 1 1 -1 -1 -1 1 ones(1,15)]';

studydir = 'C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\included_studies';
studydir2 = 'C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\zunhammer2018_studies';

% Some fixes
% ColNames(4) = [];
% Group_labels(4) = [];
% idx_vec(4) = [];
st_vec = cat(2,[1 2 3 4 4 4 5 5 5 5 6 7 8 9 10],11:25);
C = (([1 1 1 -1 -1 -1 1 1 1 -1 -1 -1 -1 1 1 -1 -1 -1 1 1 -1 -1 -1 -1 -1 1 1 -1 1 1]+1)/2)+1;
Idx = {[1], [1], [2], [4:6],[1:4],[1],[1], [1], [1], [2]}';

[BigDat_f] = harvest_canlab_funccols(studydir, Idx);
ncontrast = numel(BigDat_f);
for ii = 1:ncontrast
    % BigDat_f{ii}=  rescale(BigDat_f{ii}, 'csf_mean_var');
    BigDat_f{ii}=  rescale(BigDat_f{ii}, 'l2norm_images');
        % BigDat_f{ii}=  rescale(BigDat_f{ii}, 'prctileimages');
    BigDat_f{ii}.dat  = BigDat_f{ii}.dat.*idx_vec(ii); % Change signs of the contrasts so that all contrasts are for positive effects
end

% Zunhammer studies
Idx2 = num2cell(ones(15,1));
[BigDat_f2] = harvest_canlab_funccols(studydir2, Idx2);
ncontrast = numel(BigDat_f2);
for ii = 1:ncontrast
    % BigDat_f2{ii}=  rescale(BigDat_f2{ii}, 'csf_mean_var');
    BigDat_f2{ii}=  rescale(BigDat_f2{ii}, 'l2norm_images');
     % BigDat_f2{ii}=  rescale(BigDat_f2{ii}, 'prctileimages');
    BigDat_f2{ii}.dat  = BigDat_f2{ii}.dat.*-1; % Change signs of the contrasts so that all contrasts are for positive effects
end
BigDat_f = cat(2,BigDat_f, BigDat_f2);

% Normalize to same space:
data_cell = alignFmriDataToReference(BigDat_f, [2]);
data_cell_rn = cellfun(@(x) harmonize_zero_preserve(x),data_cell,'UniformOutput',false);
% data_cell_rn = data_cell ;

% C = C(1:16);
% st_vec = st_vec(1:16); 
[tmap_iv, pmap_iv, df] = voxelwiseLM(data_cell_rn, C, st_vec);
dsc = cellfun(@(x) x.dat',data_cell_rn,'UniformOutput',false);

%% NPS masks
ord = load('C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\Postprocessing\ratings_sort.mat');
ord = ord.ord;
Sgn_dat = apply_all_signatures(data_cell,'conditionnames',ColNames);

Group_labels = ColNames_mat.group_var;
del_idx = [4 8 10];
Group_labels(del_idx)=[];
Group_labels = Group_labels(ord);
 
cl_mat = lines(length(unique(Group_labels)));
labels = {'NPS','NPSpos','NPSneg','SIIPS','PINES','Rejection','VPS','FM_Multisens','FM_pain','Empathic_Care'};

fig = figure('Color','w','Position',[100 100 2400 1200]);
hold on
k = 0;
showlegender = false; 
showax = false;
for zz = 1:length(labels)
    k = k+1;
    ax = subplot(5,2,k);
    if zz==length(labels); showlegender = true; end
    if zz==9||zz==length(labels); showax  = true; else; showax =false; end
    eval(sprintf('cohensD = cohensD_table_wrapper(Sgn_dat.%s);',labels{zz}))

    cohensD(del_idx,:)=[];
    cohensD = cohensD(ord,:);


   % plotCohensD_byGroup(cohensD,Group_labels,'Color',cl_mat,'Sorting',false)
    % boxplot_table_sorted_by_median(T,Group_labels,'SortMode','Group')

    [figOut, order] = plotCohensD_byGroup(cohensD, Group_labels, ...
        'Axes', ax, ...         % tells function to plot into this subplot
        'Colors', cl_mat, ...        % supply full color matrix matching canonical labels
        'Sorting', false, ...
        'ShowLegend', showlegender,...
        'ShowXlabel', showax ,...
        'FigureTitle',labels{zz}); 
end

%% Clustering
% Apply parcellation
atlasname = 'canlab2024_coarse_fmriprep20_2mm';
atlas = load_atlas(atlasname);
f_correct = cell(1,30);
[contrastCells] = apply_parc_contrasiobj(data_cell_rn,atlas,f_correct);  

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
