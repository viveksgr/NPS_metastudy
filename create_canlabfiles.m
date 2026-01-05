% T is a table (nSubjects × nVars)
% Example:
% T = table(subject_ids, age, sex, RT_mean, acc_mean, 'VariableNames', {'subjid','age','sex','RT','acc'});

DAT = canlab_dataset();           % empty dataset object

% Subject-level variables
DAT.Subj_Level.names{1} = 'nodrug_vs_drug';   % cell array of names
DAT.Subj_Level.data  = WO1;               % numeric matrix (nSubjects x nVars)
% if your table contains non-numeric columns (strings/cells), convert appropriately:
% e.g. D.Subj_Level.data = [ numericCols ... ];
% and keep text vars separately in D.Subj_Level.labels (or convert to categorical numeric codes)

% optional bookkeeping
% if ismember('subjid', T.Properties.VariableNames)
    DAT.Subj_Level.subjects = wh_subjs;   % it's useful to store subject ids
% end

% Save
save('canlab_dataset_atlas_2013.mat','DAT');

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

%%%

