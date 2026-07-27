function convert_to_nii_realign_to_first()
% Convert Analyze .hdr/.img.gz pairs to NIfTI, then realign/reslice to the
% first subject in numeric subject order using SPM12.
%
% Outputs:
%   nii_unaligned/         original converted .nii files
%   nii_realign_to_first/  copied .nii inputs plus SPM r*.nii resliced files
%
% Run from MATLAB:
%   cd('C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\additional_studies_toinc\Roy2009\neg_vs_pos')
%   convert_to_nii_realign_to_first

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

if exist('spm', 'file') ~= 2
    candidateSpmDirs = { ...
        fullfile(getenv('USERPROFILE'), 'Dartmouth College Dropbox', 'Vivek Sagar', 'toolboxes', 'spm12'), ...
        fullfile(getenv('USERPROFILE'), 'Dartmouth College Dropbox', 'Vivek Sagar', 'mind', 'toolboxes', 'spm12')};
    for iSpm = 1:numel(candidateSpmDirs)
        if isfile(fullfile(candidateSpmDirs{iSpm}, 'spm.m'))
            addpath(candidateSpmDirs{iSpm});
            break
        end
    end
end

if exist('spm', 'file') ~= 2
    error(['SPM12 is not on the MATLAB path. Add SPM12 first, for example: ' ...
           'addpath(''C:\path\to\spm12'')']);
end

spm('Defaults', 'fmri');
spm_jobman('initcfg');

unalignedDir = fullfile(scriptDir, 'nii_unaligned');
realignDir = fullfile(scriptDir, 'nii_realign_to_first');
tempDir = fullfile(scriptDir, '_tmp_analyze_unzipped');

ensure_dir(unalignedDir);
ensure_dir(realignDir);
ensure_dir(tempDir);

hdrFiles = dir(fullfile(scriptDir, '*.hdr'));
if isempty(hdrFiles)
    error('No .hdr files found in %s.', scriptDir);
end

[hdrFiles, subjectNums] = sort_by_subject_number(hdrFiles);
fprintf('Found %d Analyze pairs.\n', numel(hdrFiles));
fprintf('Reference/first subject: sub%d (%s)\n', subjectNums(1), hdrFiles(1).name);

niiFiles = cell(numel(hdrFiles), 1);
for i = 1:numel(hdrFiles)
    [~, baseName] = fileparts(hdrFiles(i).name);
    hdrPath = fullfile(scriptDir, hdrFiles(i).name);
    gzPath = fullfile(scriptDir, [baseName '.img.gz']);
    if ~isfile(gzPath)
        error('Missing image file for %s: %s', hdrFiles(i).name, gzPath);
    end

    tempHdr = fullfile(tempDir, [baseName '.hdr']);
    tempImg = fullfile(tempDir, [baseName '.img']);
    outNii = fullfile(unalignedDir, [baseName '.nii']);
    niiFiles{i} = outNii;

    copyfile(hdrPath, tempHdr);
    if ~isfile(tempImg)
        gunzip(gzPath, tempDir);
    end

    if ~isfile(outNii)
        fprintf('Converting %s -> %s\n', hdrFiles(i).name, outNii);
        V = spm_vol(tempHdr);
        Y = spm_read_vols(V);
        Vout = V;
        Vout.fname = outNii;
        spm_write_vol(Vout, Y);
    else
        fprintf('Keeping existing conversion: %s\n', outNii);
    end
end

realignInputs = cell(numel(niiFiles), 1);
for i = 1:numel(niiFiles)
    [~, baseName, ext] = fileparts(niiFiles{i});
    dst = fullfile(realignDir, [baseName ext]);
    copyfile(niiFiles{i}, dst);
    realignInputs{i} = [dst ',1'];
end

matlabbatch = {};
matlabbatch{1}.spm.spatial.realign.estwrite.data = {realignInputs};
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

fprintf('Realigning/reslicing to %s...\n', realignInputs{1});
spm_jobman('run', matlabbatch);

fprintf('\nDone.\n');
fprintf('Converted NIfTIs: %s\n', unalignedDir);
fprintf('Resliced NIfTIs:  %s\n', realignDir);

end

function ensure_dir(pathName)
if ~exist(pathName, 'dir')
    mkdir(pathName);
end
end

function [sortedFiles, subjectNums] = sort_by_subject_number(files)
subjectNums = nan(numel(files), 1);
for i = 1:numel(files)
    token = regexp(files(i).name, '^sub(\d+)_', 'tokens', 'once');
    if isempty(token)
        error('Could not parse subject number from filename: %s', files(i).name);
    end
    subjectNums(i) = str2double(token{1});
end

[subjectNums, order] = sort(subjectNums);
sortedFiles = files(order);
end
