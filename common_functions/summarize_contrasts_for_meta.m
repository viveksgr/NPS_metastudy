function summarize_contrasts_for_meta()
% SUMMARIZE_CONTRASTS_FOR_META
%   1. Pull con_0001 / con_0003 (.img + .hdr, optionally .gz) from each
%      MPA2* subject in ExperienceSession and RegulateSession.
%   2. Copy into <root>/summarized_for_meta/{con0001,con0003}/{Experience,Regulate}
%      with a subject prefix in the filename.
%   3. Convert Analyze (.img/.hdr) -> NIfTI (.nii).
%   4. Compute per-subject Experience - Regulate difference for each contrast.
%   5. Coregister all per-subject difference images to a common reference
%      (first subject's con_0001 difference) so they share a voxel grid for
%      group analysis.
%
% Requires SPM12 on the MATLAB path.

    % ---- paths -----------------------------------------------------------
    root = '\\dartfs-hpc.dartmouth.edu\rc\lab\C\CANlab\labdata\projects\MPA2\Analysis\Model1\First_level';
    sessions  = {'ExperienceSession', 'RegulateSession'};
    sessTags  = {'Experience',        'Regulate'};
    contrasts = {'con_0001', 'con_0003'};
    conTags   = {'con0001',  'con0003'};

    outRoot = fullfile(root, 'summarized_for_meta');
    if ~exist(outRoot, 'dir'); mkdir(outRoot); end

    if isempty(which('spm'))
        error('SPM not found on path. addpath(genpath(...)) for SPM12 first.');
    end
    spm('defaults','FMRI');

    % ---- subject intersection across sessions ---------------------------
    subjSets = cell(1, numel(sessions));
    for s = 1:numel(sessions)
        d = dir(fullfile(root, sessions{s}, 'MPA2*'));
        d = d([d.isdir]);
        subjSets{s} = {d.name};
    end
    subjects = intersect(subjSets{1}, subjSets{2});
    fprintf('Found %d subjects present in both sessions.\n', numel(subjects));

    % ---- step 1+2+3: copy + gunzip + convert to .nii --------------------
    niiPaths = struct();   % niiPaths.(conTag).(sessTag){i}
    for c = 1:numel(contrasts)
        for s = 1:numel(sessions)
            outDir = fullfile(outRoot, conTags{c}, sessTags{s});
            if ~exist(outDir, 'dir'); mkdir(outDir); end
            niiPaths.(conTags{c}).(sessTags{s}) = cell(numel(subjects),1);
        end
    end

    for i = 1:numel(subjects)
        subj = subjects{i};
        for c = 1:numel(contrasts)
            for s = 1:numel(sessions)
                srcDir = fullfile(root, sessions{s}, subj);
                outDir = fullfile(outRoot, conTags{c}, sessTags{s});

                imgSrc = locate_pair(srcDir, contrasts{c}, '.img');
                hdrSrc = locate_pair(srcDir, contrasts{c}, '.hdr');
                if isempty(imgSrc) || isempty(hdrSrc)
                    warning('Missing %s for %s / %s', contrasts{c}, subj, sessions{s});
                    continue;
                end

                imgDst = fullfile(outDir, [subj '_' contrasts{c} '.img']);
                hdrDst = fullfile(outDir, [subj '_' contrasts{c} '.hdr']);
                copy_maybe_gunzip(imgSrc, imgDst);
                copy_maybe_gunzip(hdrSrc, hdrDst);

                % Analyze -> NIfTI via SPM
                niiDst = fullfile(outDir, [subj '_' contrasts{c} '.nii']);
                V = spm_vol(imgDst);
                Y = spm_read_vols(V);
                V.fname = niiDst;
                V.private = [];               % force re-create
                spm_write_vol(V, Y);

                niiPaths.(conTags{c}).(sessTags{s}){i} = niiDst;
            end
        end
    end

    % ---- step 4: per-subject Experience - Regulate difference -----------
    diffPaths = struct();
    for c = 1:numel(contrasts)
        diffDir = fullfile(outRoot, conTags{c}, 'diff_Exp_minus_Reg');
        if ~exist(diffDir,'dir'); mkdir(diffDir); end
        diffPaths.(conTags{c}) = cell(numel(subjects),1);

        for i = 1:numel(subjects)
            expF = niiPaths.(conTags{c}).Experience{i};
            regF = niiPaths.(conTags{c}).Regulate{i};
            if isempty(expF) || isempty(regF); continue; end

            outF = fullfile(diffDir, [subjects{i} '_' contrasts{c} '_ExpMinusReg.nii']);
            spm_imcalc({expF; regF}, outF, 'i1 - i2', ...
                       struct('interp',1,'mask',0,'dtype',spm_type('float32')));
            diffPaths.(conTags{c}){i} = outF;
        end
    end

    % ---- step 5: coregister diffs to a common reference -----------------
    % Use first available con_0001 diff as reference; coregister all other
    % diff images (both contrasts) to it. Reslice so they share grid/voxels.
    refIdx = find(~cellfun(@isempty, diffPaths.con0001), 1, 'first');
    if isempty(refIdx)
        warning('No diff images produced; skipping coregistration.');
        return;
    end
    refImg = diffPaths.con0001{refIdx};

    % Reslice (no re-estimation) all diff images onto the reference grid.
    % Con images are already in a common (normalized) space from the
    % first-level GLM — we only need a shared voxel grid, not a new
    % rigid-body transform between subjects.
    P = {refImg};
    for c = 1:numel(contrasts)
        for i = 1:numel(subjects)
            f = diffPaths.(conTags{c}){i};
            if isempty(f) || strcmp(f, refImg); continue; end
            P{end+1,1} = f; %#ok<AGROW>
        end
    end

    flags = struct('mask',0, 'mean',0, 'interp',4, 'which',1, 'wrap',[0 0 0], 'prefix','r');
    spm_reslice(P, flags);

    fprintf('Done. Reference: %s\n', refImg);
    fprintf('Resliced diffs are r*.nii alongside each source.\n');
end

% =====================================================================
function p = locate_pair(srcDir, base, ext)
    cand = { fullfile(srcDir, [base ext])      ...
             fullfile(srcDir, [base ext '.gz']) };
    p = '';
    for k = 1:numel(cand)
        if exist(cand{k}, 'file'); p = cand{k}; return; end
    end
end

function copy_maybe_gunzip(src, dst)
    if endsWith(src, '.gz')
        tmp = gunzip(src, tempdir);
        movefile(tmp{1}, dst);
    else
        copyfile(src, dst);
    end
end
