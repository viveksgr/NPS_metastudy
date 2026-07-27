function outFiles = realign_reslice_to_first(niiFiles, varargin)
%REALIGN_RESLICE_TO_FIRST Realign and reslice NIfTI files to the first image.
%
%   outFiles = realign_reslice_to_first(niiFiles)
%
% Uses SPM to estimate rigid-body alignment using the first image in
% niiFiles as the reference, then writes resliced images with the SPM
% default "r" prefix. By default, the reference image is not rewritten; the
% returned outFiles list contains the original reference image followed by
% the resliced remaining images.
%
% Inputs:
%   niiFiles  Cell array, string array, or char array of .nii filenames.
%
% Optional name/value inputs:
%   'ResliceFirst'  Default: false. If true, also writes a resliced copy of
%                   the first/reference image.
%   'Interp'        Default: 4. SPM interpolation order.
%   'Prefix'        Default: 'r'. Prefix for resliced output files.
%   'Verbose'       Default: true.
%
% Example:
%   files = {
%       'sub01_con_0001.nii'
%       'sub02_con_0001.nii'
%       'sub03_con_0001.nii'
%   };
%   outFiles = realign_reslice_to_first(files);

opts = inputParser;
opts.addRequired('niiFiles', @(x) iscell(x) || isstring(x) || ischar(x));
opts.addParameter('ResliceFirst', false, @islogical);
opts.addParameter('Interp', 4, @(x) isnumeric(x) && isscalar(x));
opts.addParameter('Prefix', 'r', @(x) ischar(x) || isstring(x));
opts.addParameter('Verbose', true, @islogical);
opts.parse(niiFiles, varargin{:});

niiFiles = normalize_file_list(opts.Results.niiFiles);
resliceFirst = opts.Results.ResliceFirst;
interp = opts.Results.Interp;
prefix = char(opts.Results.Prefix);
verbose = opts.Results.Verbose;

if numel(niiFiles) < 2
    error('Provide at least two .nii files to realign/reslice.');
end

if exist('spm', 'file') ~= 2
    error('SPM is not on the MATLAB path. Add SPM before running this function.');
end

for iFile = 1:numel(niiFiles)
    thisFile = char(niiFiles(iFile));
    [~, ~, ext] = fileparts(strip_volume_index(thisFile));

    if ~strcmpi(ext, '.nii')
        error('Expected a .nii file, got: %s', thisFile);
    end

    if ~isfile(strip_volume_index(thisFile))
        error('File does not exist: %s', thisFile);
    end
end

spm('defaults', 'fmri');

spmFiles = cellstr(add_volume_index(niiFiles));
P = char(spmFiles);

realignFlags = struct( ...
    'quality', 0.9, ...
    'sep', 4, ...
    'fwhm', 5, ...
    'rtm', 0, ...
    'interp', 2, ...
    'wrap', [0 0 0], ...
    'weight', '');

resliceFlags = struct( ...
    'mask', 1, ...
    'mean', 0, ...
    'interp', interp, ...
    'which', 1 + double(resliceFirst), ...
    'wrap', [0 0 0], ...
    'prefix', prefix);

if verbose
    fprintf('Realigning %d images to reference: %s\n', numel(niiFiles), niiFiles(1));
end

spm_realign(P, realignFlags);
spm_reslice(P, resliceFlags);

outFiles = strings(numel(niiFiles), 1);
for iFile = 1:numel(niiFiles)
    if iFile == 1 && ~resliceFirst
        outFiles(iFile) = niiFiles(iFile);
        continue;
    end

    [folderPath, fileName, ext] = fileparts(char(niiFiles(iFile)));
    outFiles(iFile) = string(fullfile(folderPath, [prefix fileName ext]));
end
end

function files = normalize_file_list(filesIn)
if ischar(filesIn)
    files = string(cellstr(filesIn));
elseif iscell(filesIn)
    files = string(filesIn(:));
else
    files = string(filesIn(:));
end

files = strtrim(files);
files(files == "") = [];

if isempty(files)
    error('No input files were provided.');
end
end

function files = add_volume_index(files)
filesIn = files;
files = strings(size(filesIn));
for iFile = 1:numel(filesIn)
    fileName = char(filesIn(iFile));
    if isempty(regexp(fileName, ',\d+$', 'once'))
        files(iFile) = string([fileName ',1']);
    else
        files(iFile) = string(fileName);
    end
end
end

function fileName = strip_volume_index(fileName)
fileName = regexprep(char(fileName), ',\d+$', '');
end
