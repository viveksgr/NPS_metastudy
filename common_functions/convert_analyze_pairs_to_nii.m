function convert_analyze_pairs_to_nii(inputDir, outputDir, overwriteExisting)
% Convert Analyze .hdr/.img.gz image pairs in a folder to single-file .nii.
%
% Usage:
%   cd('C:\Users\sgrvi\Dartmouth College Dropbox\Vivek Sagar\Sagar_2025_Pain_Intervention_Meta_Analysis_PIMA\Data\subjectlevel\additional_studies_toinc\Roy2009\neg_vs_pos')
%   convert_analyze_pairs_to_nii
%
% Optional:
%   convert_analyze_pairs_to_nii(inputDir)
%   convert_analyze_pairs_to_nii(inputDir, outputDir)
%   convert_analyze_pairs_to_nii(inputDir, outputDir, true)  % overwrite

if nargin < 1 || isempty(inputDir)
    inputDir = fileparts(mfilename('fullpath'));
end

if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(inputDir, 'nii');
end

if nargin < 3 || isempty(overwriteExisting)
    overwriteExisting = false;
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

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

tempDir = fullfile(inputDir, '_tmp_analyze_unzipped');
if ~exist(tempDir, 'dir')
    mkdir(tempDir);
end

hdrFiles = dir(fullfile(inputDir, '*.hdr'));
if isempty(hdrFiles)
    error('No .hdr files found in %s.', inputDir);
end

hdrFiles = sort_by_subject_number(hdrFiles);
fprintf('Found %d Analyze header files in %s.\n', numel(hdrFiles), inputDir);

for iFile = 1:numel(hdrFiles)
    [~, baseName] = fileparts(hdrFiles(iFile).name);
    hdrPath = fullfile(inputDir, [baseName '.hdr']);
    gzPath = fullfile(inputDir, [baseName '.img.gz']);
    tempHdrPath = fullfile(tempDir, [baseName '.hdr']);
    tempImgPath = fullfile(tempDir, [baseName '.img']);
    niiPath = fullfile(outputDir, [baseName '.nii']);

    if ~isfile(gzPath)
        error('Missing matching .img.gz for %s.', hdrPath);
    end

    if isfile(niiPath) && ~overwriteExisting
        fprintf('Skipping existing file: %s\n', niiPath);
        continue
    end

    copyfile(hdrPath, tempHdrPath);
    if ~isfile(tempImgPath) || overwriteExisting
        if isfile(tempImgPath)
            delete(tempImgPath);
        end
        gunzip(gzPath, tempDir);
    end

    fprintf('Converting %s -> %s\n', hdrFiles(iFile).name, niiPath);
    V = spm_vol(tempHdrPath);
    Y = spm_read_vols(V);
    Vout = V;
    Vout.fname = niiPath;
    spm_write_vol(Vout, Y);
end

fprintf('Done. NIfTI files are in: %s\n', outputDir);

end

function sortedFiles = sort_by_subject_number(files)
subjectNums = nan(numel(files), 1);
for iFile = 1:numel(files)
    token = regexp(files(iFile).name, '^sub(\d+)_', 'tokens', 'once');
    if isempty(token)
        subjectNums(iFile) = iFile;
    else
        subjectNums(iFile) = str2double(token{1});
    end
end

[~, order] = sort(subjectNums);
sortedFiles = files(order);
end
