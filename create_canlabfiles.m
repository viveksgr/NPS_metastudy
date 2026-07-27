% T is a table (nSubjects × nVars)
% Example:
% T = table(subject_ids, age, sex, RT_mean, acc_mean, 'VariableNames', {'subjid','age','sex','RT','acc'});

DAT = canlab_dataset();           % empty dataset object

% Add ratings manually from excel
subj_dat = [];
subj_nam = [];

% Subject-level variables
DAT.Subj_Level.names{1} = 'nature-indoor';   % cell array of names
DAT.Subj_Level.data  = subj_dat;               % numeric matrix (nSubjects x nVars)
% if your table contains non-numeric columns (strings/cells), convert appropriately:
% e.g. D.Subj_Level.data = [ numericCols ... ];
% and keep text vars separately in D.Subj_Level.labels (or convert to categorical numeric codes)
DAT.Subj_Level.subjects = subj_nam;   % it's useful to store subject ids

mkdir('canlab_dataset_objects_for_pain_ratings')
% Save
save('canlab_dataset_objects_for_pain_ratings\canlab_dataset','DAT');

% Identify nums
% Atlas 2013
fx = @(x) sscanf(x(1:8),'%*[^0-9]%d');
cnames_17 = cellfun(fx,DAT_set.Subj_Level.id);
% fx = @(x) sscanf(x,'%*[^0-9]%d');
cnames_19 = cellfun(fx,wh_subjs);

nodrg = [nanmean(HH); nanmean(HO)];
drg = [nanmean(WH); nanmean(WO)];
nodrg = nanmean(nodrg);
drg = nanmean(drg);
diff = nodrg'-drg';
find(~lia)
diff([11 19]) = [];

HH1 = nanmean(HH);
HH1([11 19])=[];
HH1 = HH1';

HO1 = nanmean(HO);
HO1([11 19])=[];
HO1 = HO1';

WH1 = nanmean(WH);
WH1([11 19])=[];
WH1 = WH1';

WO1 = nanmean(WO);
WO1([11 19])=[];
WO1 = WO1';

wo = -(WH1+WO1-HH1-HO1);
wh_subjs([11 19])=[];

%% reslice

outFiles = realign_reslice_to_first(files, 'ResliceFirst', true);


flags = struct('mask',0,'mean',0,'interp',1,'which',1,'wrap',[0 0 0],'prefix','r');

confiles = dir(fullfile(pwd,'sub*.nii')); % adjust
% % build fullpaths cell

newname = strcat('r',confiles(1).name);
newpath = fullfile(confiles(1).folder,newname);
confiles = fullfile({confiles.folder}, {confiles.name});
spm_reslice(confiles, flags);

copyfile(confiles{1},newpath) % Copy the first file w/o modification 

%% realign

confiles = dir(fullfile(pwd,'sub*.nii')); % adjust
confiles = fullfile({confiles.folder}, {confiles.name});
outFiles = realign_reslice_to_first(confiles, 'ResliceFirst', true);
% % build fullpaths cell

%% contrast obj

% load('data_frame.mat')
% gather all subject contrast files (example pattern)
confiles = dir(fullfile(pwd,'r*007.nii')); % adjust
% % build fullpaths cell
confiles = fullfile({confiles.folder}, {confiles.name});
% create fmri_data object (uses CanlabCore fmri_data constructor)
fmriDat = fmri_data(confiles);    % will apply default brain mask (or specify mask)

% c = load_image_set('nps');
% fmriDat2 = resample_space(fmriDat,c);

% Sgn_dat2 = apply_nps(fmriDat);
Sgn_dat2 = apply_siips(fmriDat);
nps = Sgn_dat2{1};
npsp = nps/std(nps,0);
mean(npsp)
1.96*std(npsp)/sqrt(length(npsp))

fmriDat.dat = fmriDat.dat;
fmriDat.source_notes = 'predict-uncontrol';
DATA_OBJ_CON{1} = fmriDat;
mkdir("fmri_data_objects_for_contrasts")
save("fmri_data_objects_for_contrasts\contrast_data_objects.mat","DATA_OBJ_CON")

%%
Sgn_dat2 = apply_nps(DATA_OBJ_CON{1});
nps = Sgn_dat2{1};
mean(nps)/std(nps,0)

%% Compile single trials
addpath('\\dartfs-hpc.dartmouth.edu\rc\lab\C\CANlab\labdata\projects\canlab_single_trials_for_git_repo')

% Paingen
N_heat = image_obj.metadata_table.heat;
N_plac = image_obj.metadata_table.prodicaine;
N_plac(~logical(N_heat))=-1;
N_sub = image_obj.metadata_table.subject_id;
N_rat = image_obj.metadata_table.intRating;

% ILCP0ctrl
name_c = 'Ctrl-NoCtrl';
image_obj = ilcp;
N_plac = image_obj.metadata_table.ctrl;
N_sub = image_obj.metadata_table.subject_id;
N_rat = image_obj.metadata_table.rating;

% placebo
name_c = 'placebo';
image_obj = tmpData;
N_plac = image_obj.metadata_table.placebo;
N_sub = image_obj.metadata_table.subject_id;
N_rat = image_obj.metadata_table.rating;

% Scebl/ie
name_c = 'reg_persp';
image_obj = data;
cue = image_obj.metadata_table.reg;
N_plac = 1-((image_obj.metadata_table.cue)+1/2); % 1 (-1) = low, 0 (1) = high
N_sub = image_obj.metadata_table.subject_id;
N_rat = image_obj.metadata_table.rating;

[beh_delta, subj_ids, beh_stats] = subject_condition_delta(N_rat, N_sub, N_plac);
[fmri_delta, subj_ids, fmri_stats] = subject_condition_delta_matrix(image_obj.dat,  N_sub, N_plac);

image_obj.dat = fmri_delta;
fmriDat = image_obj;
fmriDat.source_notes = name_c;
[nvox,nsub] = size(fmriDat.dat);
fmriDat.removed_images = false(nsub,1);
fmriDat.removed_voxels = false(nvox,1);
DATA_OBJ_CON{1} = fmriDat;
mkdir("fmri_data_objects_for_contrasts")
save("fmri_data_objects_for_contrasts\contrast_data_objects.mat","DATA_OBJ_CON")

DAT = canlab_dataset();           % empty dataset object
% Subject-level variables
DAT.Subj_Level.names{1} = name_c;   % cell array of names
DAT.Subj_Level.data  = beh_delta';               % numeric matrix (nSubjects x nVars)
% if your table contains non-numeric columns (strings/cells), convert appropriately:
% e.g. D.Subj_Level.data = [ numericCols ... ];
% and keep text vars separately in D.Subj_Level.labels (or convert to categorical numeric codes)
DAT.Subj_Level.subjects = subj_ids;   % it's useful to store subject ids

mkdir("canlab_dataset_objects_for_pain_ratings")
save("canlab_dataset_objects_for_pain_ratings\canlab_dataset.mat","DAT")

% concat_fmriobj(folder1, folder2)