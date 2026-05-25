% summarize_ilcp_contrasts.m
% 1) Pull con_0030 and con_0031 from each ilcp* subject folder
% 2) Convert each Analyze (.img/.hdr) pair into a single .nii in
%    summarized_for_meta/con_0030/ or summarized_for_meta/con_0031/
%    with the subject folder name as prefix
% 3) Reslice all images within each contrast folder to a common voxel grid

addpath('C:\Work\Toolboxes_general\spm12');
spm('defaults','FMRI');

ilcpDir = ['\\dartfs-hpc.dartmouth.edu\rc\lab\C\CANlab\labdata\archive_pub\2017_Woo_SIIPS_Schmidt_ILCP\Imaging\Analyses\first_level\model5'];

outDir = fullfile(ilcpDir, 'summarized_for_meta');
if ~exist(outDir,'dir'); mkdir(outDir); end

contrasts = {'con_0040','con_0041'};

% Create per-contrast subfolders
for c = 1:numel(contrasts)
    d = fullfile(outDir, contrasts{c});
    if ~exist(d,'dir'); mkdir(d); end
end

subjFolders = dir(fullfile(ilcpDir,'ilcp*'));
subjFolders = subjFolders([subjFolders.isdir]);

%% Step 1+2: gunzip .img.gz, read Analyze pair, write single .nii
copiedImgs = struct();
for c = 1:numel(contrasts); copiedImgs.(contrasts{c}) = {}; end

for s = 1:numel(subjFolders)
    sName = subjFolders(s).name;
    sDir  = fullfile(ilcpDir, sName);

    for c = 1:numel(contrasts)
        con = contrasts{c};
        srcImgGz = fullfile(sDir, [con '.img.gz']);
        srcHdr   = fullfile(sDir, [con '.hdr']);

        if ~exist(srcImgGz,'file') || ~exist(srcHdr,'file')
            warning('Missing %s for %s — skipping', con, sName);
            continue;
        end

        % Gunzip alongside header in a temp staging dir
        stage = fullfile(tempdir, ['ilcp_stage_' sName '_' con]);
        if ~exist(stage,'dir'); mkdir(stage); end
        gunzip(srcImgGz, stage);
        copyfile(srcHdr, fullfile(stage, [con '.hdr']));

        % Read Analyze pair, write single .nii into per-contrast subfolder
        V = spm_vol(fullfile(stage, [con '.img']));
        Y = spm_read_vols(V);

        dstNii = fullfile(outDir, con, sprintf('%s_%s.nii', sName, con));
        Vout       = V;
        Vout.fname = dstNii;
        spm_write_vol(Vout, Y);

        rmdir(stage,'s');

        copiedImgs.(con){end+1,1} = dstNii;
    end
end

fprintf('Wrote %d subjects worth of .nii contrasts to %s\n', numel(subjFolders), outDir);

%% Step 3: reslice within each contrast folder to a common grid
flags = struct('mask',0,'mean',0,'interp',1,'which',1,'wrap',[0 0 0],'prefix','r');
for c = 1:numel(contrasts)
    imgs = copiedImgs.(contrasts{c});
    if numel(imgs) < 2; continue; end
    spm_reslice(imgs, flags);
    fprintf('Resliced %d images for %s\n', numel(imgs), contrasts{c});
end

fprintf('Done. Output in %s\n', outDir);
