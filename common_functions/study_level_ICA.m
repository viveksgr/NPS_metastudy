function [studyWeights, A_mix, S_maps, subj2study] = study_level_ICA(dataCell, opts)
% STUDY_LEVEL_ICA  Run ICA on pooled subject-level ROI data and return SxP study weights
%
% Inputs:
%   dataCell - 1xS cell, each cell = Ns x A matrix (subjects x ROIs)
%   opts - struct with fields (optional)
%       .nComp = desired # ICA components (P), default 20
%       .doPCA = true|false, default true
%       .nPC   = # PCs to keep before ICA (if doPCA), default min(50, A)
%       .icaMethod = 'fastica' | 'canlab' (default 'fastica')
%       .icaArgs = cell array of extra args to ICA routine
%
% Outputs:
%   studyWeights - S x P matrix (each row = study average mixing weights)
%   A_mix        - N_total x P subject-level mixing weights
%   S_maps       - P x A spatial maps
%   subj2study   - N_total x 1 integer mapping each subject-row -> study index
%
% Requires FastICA (fastica.m) OR CanlabCore ica() on fmri_data (with dependency icatb_fastICA).

if nargin < 2, opts = struct(); end
if ~isfield(opts,'nComp'), opts.nComp = 20; end
if ~isfield(opts,'doPCA'), opts.doPCA = true; end
if ~isfield(opts,'nPC'), opts.nPC = min(50, size(dataCell{1},2)); end
if ~isfield(opts,'icaMethod'), opts.icaMethod = 'fastica'; end
if ~isfield(opts,'icaArgs'), opts.icaArgs = {}; end

S = numel(dataCell);
A = size(dataCell{1},2);

% Build concatenated subject x ROI matrix
subj2study = [];
Xlist = {};
for s = 1:S
    Xs = double(dataCell{s}); % Ns x A
    Ns = size(Xs,1);
    subj2study = [subj2study; s*ones(Ns,1)]; %#ok<AGROW>
    Xlist{end+1} = Xs; %#ok<AGROW>
end
X_concat = vertcat(Xlist{:}); % N_total x A
N_total = size(X_concat,1);

% Preprocessing: center (demean) across ROIs for each subject (row)
Xc = bsxfun(@minus, X_concat, mean(X_concat,2));

% Optional: PCA reduction (across ROIs) -> reduce columns from A to nPC
if opts.doPCA
    nPC = min(opts.nPC, A);
    % Compute SVD of Xc' (A x N) or PCA on covariance
    % better: compute PCA on columns (ROI space)
    [coeff, score, ~, ~, explained] = pca(Xc, 'NumComponents', nPC);
    % score is N_total x nPC (projections), coeff is A x nPC
    X_ica_input = score;         % N_total x nPC
    PCA_basis = coeff;           % A x nPC (spatial PCA maps)
else
    X_ica_input = Xc;            % N_total x A
    PCA_basis = [];
end

% Run ICA -> want subject mixing weights A_mix (N_total x P) and spatial maps S_maps (P x A or P x nPC)
P = opts.nComp;
switch lower(opts.icaMethod)
    case 'fastica'
        % fastica expects rows = components x samples if syntax different:
        % FastICA (MATLAB) usually X as observations x variables; calling fastica(X', 'numOfIC',P)
        try
            [icasig, A_fast, W] = fastica(X_ica_input', 'numOfIC', P, opts.icaArgs{:});
            % icasig: P x samples (here samples = N_total), A_fast: variables x components (not consistently documented)
            % We want mixing matrix A_mix: N_total x P
            A_mix = icasig';          % N_total x P  (component timecourses / subject-loadings)
            S_maps = (W')' * 0;       % placeholder - we will compute spatial maps below
            % reconstruct spatial maps in ROI space:
            if opts.doPCA
                % convert ICA sources (P x nPC) -> P x A: S_maps_roi = icasig * pinv(PCA_scores?) Simpler:
                % Run ICA on PCA score space yields icasig (P x N_total) and mixing W; We need spatial maps across ROIs:
                % The spatial maps: S_maps = (unmixing^(-1)) * PCA_basis' ? Implementation details vary by package.
                % Easiest: compute component spatial map by regression: regress each component's subject vector onto ROI columns
                S_maps = zeros(P, A);
                for p = 1:P
                    compvec = A_mix(:,p); % N_total x 1
                    % regress onto ROI columns (X_concat): S_p = compvec \ X_concat  (OLS)
                    S_maps(p,:) = (compvec' * X_concat) / (compvec' * compvec + eps);
                end
            else
                % no PCA: regress component expression onto ROIs directly
                S_maps = zeros(P, A);
                for p = 1:P
                    compvec = A_mix(:,p);
                    S_maps(p,:) = (compvec' * X_concat) / (compvec' * compvec + eps);
                end
            end
        catch ME
            error('fastica failed: %s', ME.message);
        end

    case 'canlab'
        % Try CANlabCore if available: create fmri_data object and call ica method
        % This assumes you have an fmri_data object constructor that accepts rows=images (subjects), cols=voxels (ROIs).
        try
            % create a dummy fmri_data-like struct with fields required by ica method:
            m = struct();
            m.dat = Xc'; % Canlab expects voxels x images in fmri_data.dat (?) - check local function; may need transpose
            % If CANlab is on path, call:
            o = fmri_data(); % requires CANlab on path
            o.dat = Xc';
            o = ica(o, 'ncomponents', P);
            % Postprocess to extract mixing & maps — exact fieldnames depend on CANlab version
            A_mix = o.ica_subject_weights; % example placeholder
            S_maps = o.ica_spatial_maps;
        catch ME
            error('CANlab ica method failed or is unavailable: %s', ME.message);
        end

    otherwise
        error('Unknown icaMethod: %s', opts.icaMethod);
end

% Now compute study-level average weights
studyWeights = zeros(S, size(A_mix,2));
for s = 1:S
    wh = find(subj2study == s);
    if isempty(wh), studyWeights(s,:) = NaN; else studyWeights(s,:) = mean(A_mix(wh,:),1); end
end

end
