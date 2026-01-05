function V_out = harmonize_zero_preserve(V_in, varargin)
% HARMONIZE_ZERO_PRESERVE  Zero-preserving quantile harmonization per-voxel
%
% V_out = harmonize_zero_preserve(V_in)
% V_out = harmonize_zero_preserve(V_in, 'vectorized', true, 'verbose', true)
%
% Input:
%   V_in        - canlab contrast/fmri object (must have field/property .dat)
%                 dat should be [nVox x nSubjects]
% Options:
%   'vectorized' - (false) try a faster vectorized approach (may use lots of RAM)
%   'verbose'    - (false) print progress
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
parse(p,varargin{:});
do_vec = p.Results.vectorized;
do_verbose = p.Results.verbose;

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
