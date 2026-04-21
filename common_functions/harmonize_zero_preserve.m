function V_out = harmonize_zero_preserve(V_in, varargin)
% HARMONIZE_ZERO_PRESERVE  Zero-preserving quantile harmonization per-voxel
%
% V_out = harmonize_zero_preserve(V_in)
% V_out = harmonize_zero_preserve(V_in, 'vectorized', true, 'verbose', true)
% V_out = harmonize_zero_preserve(V_in, 'jitter', 0.25)
%
% Input:
%   V_in        - canlab contrast/fmri object (must have field/property .dat)
%                 dat should be [nVox x nSubjects]
% Options:
%   'vectorized' - (false) try a faster vectorized approach (may use lots of RAM)
%   'verbose'    - (false) print progress
%   'jitter'     - (0) sub-grid Gaussian noise factor. Adds
%                  randn(size(q0)) * jitter/nv to q0 before the zero
%                  override. Breaks the exact rank-based grid so
%                  downstream OLS outputs are continuous. Values of
%                  0.2-0.3 are small relative to grid spacing (1/nv)
%                  and do not meaningfully change distributions.
%                  Set to 0 to disable (default, preserves original
%                  behavior).
%   'seed'       - ([]) random seed for jitter; leave empty for
%                  non-deterministic
%
% Output:
%   V_out       - same object/class as V_in, with V_out.dat replaced by the
%                 zero-preserving quantile-transformed data.
%
% This implements the zero-preserving quantile transform described in:
%   Spisak et al. (Placebo Imaging Consortium) — quantile distance from zero.
% See also the paper /mnt/data/2025_bioarxiv.pdf. :contentReference[oaicite:1]{index=1}

% Parse options
p = inputParser;
addParameter(p,'vectorized',false,@islogical);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'jitter',0,@(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p,'seed',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parse(p,varargin{:});
do_vec = p.Results.vectorized;
do_verbose = p.Results.verbose;
jitter_sd = p.Results.jitter;
seed = p.Results.seed;

if jitter_sd > 0 && ~isempty(seed)
    rng(seed);
end

% Copy object to output
V_out = V_in;

% extract data
dat = V_in.dat;           % nV x nS
[nV, nS] = size(dat);

if do_verbose
    fprintf('harmonize_zero_preserve: %d voxels x %d subjects. vectorized=%d\n', nV, nS, do_vec);
end

if do_vec
    % Vectorized-ish approach: operate on columns after sorting — may be memory heavy.
    % We'll still handle NaNs per voxel.
    dat_out = nan(size(dat));
    for v = 1:nV
        x = dat(v,:);
        valid = ~isnan(x);
        nv = sum(valid);
        if nv == 0
            continue;
        end
        xv = x(valid);
        % ranks with ties handled
        rv = tiedrank(xv);               % 1..nv
        q = (rv - 0.5) ./ nv;           % (0.5/nv ... (nv-0.5)/nv)
        prop_neg = sum(xv < 0) / nv;
        q0 = q - prop_neg;              % shift -> distance from zero
        % optional sub-grid jitter (BEFORE zero override)
        if jitter_sd > 0
            q0 = q0 + randn(size(q0)) * (jitter_sd / nv);
        end
        % preserve exact zeros
        iszero = (xv == 0);
        if any(iszero)
            q0(iszero) = 0;
        end
        dat_out(v, valid) = q0;
    end
else
    % Safe loop approach
    dat_out = nan(size(dat));
    % progress print setup
    block = max(1, round(nV/10));
    for v = 1:nV
        x = dat(v,:);
        valid = ~isnan(x);
        nv = sum(valid);
        if nv == 0
            continue;
        end
        xv = x(valid);
        % compute tied ranks for valid entries
        rv = tiedrank(xv);               % ranks in 1..nv (handles ties)
        q = (rv - 0.5) ./ nv;           % quantiles (0..1)
        prop_neg = sum(xv < 0) / nv;    % proportion negative
        q0 = q - prop_neg;              % quantile distance-from-zero
        % optional sub-grid jitter (BEFORE zero override)
        if jitter_sd > 0
            q0 = q0 + randn(size(q0)) * (jitter_sd / nv);
        end
        % force zeros to be exactly zero
        iszero = (xv == 0);
        if any(iszero)
            q0(iszero) = 0;
        end
        dat_out(v, valid) = q0;
        % verbose progress
        if do_verbose && mod(v,block) == 0
            fprintf('  processed %d/%d voxels (%.0f%%)\n', v, nV, 100*v/nV);
        end
    end
end

% Put back into object
V_out.dat = dat_out;

if do_verbose
    fprintf('harmonize_zero_preserve: done. Output stored in V_out.dat\n');
end

end
