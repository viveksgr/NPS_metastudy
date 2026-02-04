function [V_out, Bmap, Xused] = regress_out_wm_csf(V_in, nuisance, varargin)
% REGRESS_OUT_WM_CSF  Regress nuisance covariates (WM/CSF comps) from contrast object
%
% Usage:
%   V_out = regress_out_wm_csf(V_in)                        % uses extract_wm_csf_comps
%   V_out = regress_out_wm_csf(V_in, nuisanceMatrix)       % nuisanceMatrix: nSub x p
%   [V_out, Bmap, Xused] = regress_out_wm_csf(...,'robust',true,'addIntercept',true)
%
% Inputs:
%   V_in     - canlab contrast/fmri object (expects V_in.dat is nVox x nSub)
%   nuisance - optional: nSub x p matrix OR function handle that returns nSub x p
%              if empty or omitted, attempts to call extract_wm_csf_comps(V_in) or (V_in.dat)
%
% Name/value options:
%   'addIntercept' (true)   - prepend a column of ones to X
%   'robust'       (false)  - use robustfit per-voxel (slow)
%   'verbose'      (true)
%
% Outputs:
%   V_out   - same class as V_in, with V_out.dat = residuals (nVox x nSub)
%   Bmap    - p×nV matrix of betas (excluding intercept if addIntercept true then p includes it)
%   Xused   - the design matrix used (nSub x p_used)
%
% Notes:
%  - Handles NaNs per voxel: for voxels with missing observations we do per-voxel LS
%  - For voxels with full data (no NaNs) we solve in a vectorized way: B = X \ Y
%  - Robust option uses MATLAB robustfit (intercept handling is automatic there)
%
% Example:
%   Vclean = regress_out_wm_csf(V, [], 'addIntercept', true);

% ----------------- parse options -----------------
p = inputParser;
addParameter(p,'addIntercept',true,@islogical);
addParameter(p,'robust',false,@islogical);
addParameter(p,'verbose',true,@islogical);
parse(p,varargin{:});
addIntercept = p.Results.addIntercept;
doRobust = p.Results.robust;
verbose = p.Results.verbose;

% ----------------- load data and check dims -----------------
dat = V_in.dat; % expect nVox x nSub
if size(dat,1) < size(dat,2)
    % still OK: we require nVox x nSub -- if user passed transposed object, try to detect
    % but don't force transpose blindly. We'll check typical orientation.
    % Heuristic: if more subjects than voxels then probably transposed -> transpose
    if size(dat,2) > size(dat,1)
        % leave as-is (nVox < nSub) — user probably has fewer voxels (parcellation)
        % but we require rows = voxels. So ensure: rows = voxels.
    end
end

[nV, nS] = size(dat);
if nS < 2
    error('Input must have at least two subjects (columns). Found %d.', nS);
end

% ----------------- determine nuisance matrix -----------------
X = [];
if nargin < 2 || isempty(nuisance)
    if verbose, fprintf('No nuisance provided — calling extract_wm_csf_comps...\n'); end
    % try function call: first try as function taking object, else dat
    try
        X = extract_wm_csf_comps(V_in); % expected nSub x p
    catch
        try
            X = extract_wm_csf_comps(dat); % some versions expect raw matrix
        catch ME
            error('Could not compute nuisance matrix automatically. Provide nuisance matrix as second arg. Last error: %s', ME.message);
        end
    end
elseif isa(nuisance,'function_handle')
    % call user-supplied function (try both signatures)
    try
        X = nuisance(V_in);
    catch
        X = nuisance(dat);
    end
elseif isnumeric(nuisance)
    X = nuisance;
else
    error('nuisance must be empty, a function handle, or numeric matrix nSub x p.');
end

% Validate X size
if size(X,1) ~= nS
    error('Nuisance matrix must have rows equal to number of subjects (%d). Provided has %d rows.', nS, size(X,1));
end

% Add intercept if requested
if addIntercept
    Xused = [ones(nS,1), X];  % nSub x p_used
else
    Xused = X;
end

p_used = size(Xused,2);

% center nuisance columns? (user can do this before passing if desired)
if verbose
    fprintf('Regression design: %d predictors (including intercept=%d). Solving for %d voxels x %d subjects.\n', p_used, addIntercept, nV, nS);
end

% ----------------- pre-allocate outputs -----------------
Y = dat';                % nS x nV  (subjects x voxels)
V_out = V_in;
V_out.dat = nan(size(dat)); % fill later
Bmap = nan(p_used, nV);     % store betas for each voxel

% Identify voxels with no NaNs across subjects (fast vectorized solve)
fullCols = all(~isnan(Y), 1);   % 1 x nV logical
cols_full_idx = find(fullCols);
cols_partial_idx = find(~fullCols);

% Vectorized solve for full columns
if ~isempty(cols_full_idx)
    Yfull = Y(:, cols_full_idx);    % nS x nFull
    % Solve least squares for many voxels at once:
    % B = X \ Y  -> p_used x nFull
    Bfull = Xused \ Yfull;          % uses MATLAB backslash (handles rank-def)
    % residuals
    Efull = Yfull - Xused * Bfull;
    % assign
    V_out.dat(cols_full_idx,:) = Efull';   % nVox x nS (transposed)
    Bmap(:, cols_full_idx) = Bfull;
    if verbose
        fprintf('Vectorized solve for %d voxels (no NaNs) completed.\n', numel(cols_full_idx));
    end
end

% Per-voxel solve for columns with NaNs or partially missing data
if ~isempty(cols_partial_idx)
    if verbose, fprintf('Performing per-voxel solves for %d voxels with missing data.\n', numel(cols_partial_idx)); end
    for ii = 1:numel(cols_partial_idx)
        v = cols_partial_idx(ii);
        yv = Y(:,v);                    % nS x 1 (has NaNs)
        valid = ~isnan(yv);
        nvalid = sum(valid);
        if nvalid == 0
            % leave as NaN
            Bmap(:,v) = NaN;
            V_out.dat(v, :) = NaN;
            continue;
        end
        Xv = Xused(valid, :);  % nvalid x p_used
        yv_valid = yv(valid);
        % If robust requested, use robustfit (note robustfit includes intercept)
        if doRobust
            if addIntercept
                % robustfit expects X without intercept and adds it internally
                [b_r, ~] = robustfit(Xv(:,2:end), yv_valid); % includes intercept
                % format b_r into p_used vector: intercept first then others
                b = b_r(:);
            else
                % if no intercept requested, call robustfit with no intercept: use 'robustfit' but remove intercept
                [btmp, ~] = robustfit(Xv, yv_valid);
                % btmp includes intercept; drop it if addIntercept==false (we asked no intercept - difficult with robustfit)
                % fallback: use lsq for now
                b = (Xv \ yv_valid);
            end
        else
            % OLS via backslash on reduced rows
            b = Xv \ yv_valid;  % p_used x 1
        end
        % residuals for valid rows
        e_valid = yv_valid - Xv * b;
        % fill into output residual vector
        yres = nan(nS,1);
        yres(valid) = e_valid;
        V_out.dat(v, :) = yres';
        Bmap(:, v) = b;
    end
end

% done
if verbose
    fprintf('Finished regression. Output residuals stored in V_out.dat (nVox x nSub).\n');
end

end
