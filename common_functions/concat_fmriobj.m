function concat_fmriobj(folder1, folder2)
l1 = load(fullfile(folder1,"canlab_dataset_objects_for_pain_ratings/canlab_dataset.mat"));
DAT = l1.DAT;
l2 = load(fullfile(folder2,"canlab_dataset_objects_for_pain_ratings/canlab_dataset.mat"));
DAT2 = l2.DAT;
DAT.Subj_Level.names(2)= DAT2.Subj_Level.names(1);
DAT.Subj_Level.data(:,2)= DAT2.Subj_Level.data(:,1);

save(fullfile(folder1,"canlab_dataset_objects_for_pain_ratings/canlab_dataset.mat"),"DAT")

l2 = load(fullfile(folder2,"fmri_data_objects_for_contrasts/contrast_data_objects.mat"));
l2d2 = l2.DATA_OBJ_CON;
load(fullfile(folder1,"fmri_data_objects_for_contrasts/contrast_data_objects.mat"));
DATA_OBJ_CON(2)=l2d2(1);
save(fullfile(folder1,"fmri_data_objects_for_contrasts\contrast_data_objects.mat"),"DATA_OBJ_CON")
