% Check NPS
cd('C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\zunhammer2018_studies\Bingel_et_al_2006\raw_files2')

confiles = dir(fullfile(pwd,'r*.nii.gz')); % adjust
% % build fullpaths cell
confiles = fullfile({confiles.folder}, {confiles.name});
% create fmri_data object (uses CanlabCore fmri_data constructor)
fmriDat = fmri_data(confiles);    % will apply default brain mask (or specify mask)

nps_scores_cell = apply_nps(fmriDat);
nps_scores = nps_scores_cell{1};
nps_scores_d = nps_scores/std(nps_scores,0);
mean_effect = mean(nps_scores_d)
range = 1.96*std(nps_scores_d)/sqrt(length(nps_scores_d))