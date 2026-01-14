function out = compute_effect_size(x, varargin)
% COMPUTE_EFFECT_SIZE  compute Cohen's d, Hedges' g, SE and CIs
%
% Usage:
%   out = compute_effect_size(x, 'type','one')         % one-sample (x is vector)
%   out = compute_effect_size([pre,post], 'type','paired') % paired: x is Nx2 matrix [pre post]
%   out = compute_effect_size({x1,x2}, 'type','independent') % independent groups: cell array
%
% Options (name,value):
%   'type'   : 'one' (default) | 'paired' | 'independent'
%   'nBoot'  : 2000  bootstrap samples for CI (default 2000)
%   'alpha'  : 0.05  two-sided
%   'verbose': false
%
% Output struct with fields:
%   .d, .g, .SE_d, .CI_analytic, .CI_bootstrap, .df, .n

% parse
p = inputParser;
addParameter(p,'type','one',@(s) ismember(s,{'one','paired','independent'}));
addParameter(p,'nBoot',2000,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'alpha',0.05,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'verbose',false,@islogical);
parse(p,varargin{:});
typ = p.Results.type; nBoot = p.Results.nBoot; alpha = p.Results.alpha; verbose = p.Results.verbose;

rng('default');

switch typ
    case 'one'
        xvec = x(:);
        n = numel(xvec);
        m = mean(xvec,'omitnan');
        s = std(xvec,0,'omitnan');     % sample sd (n-1)
        d = m / s;
        df = n-1;
        % analytic SE
        SE = sqrt( (1/n) + (d.^2 ./ (2*n)) );
        % bootstrap: resample subjects
        boot_d = nan(nBoot,1);
        for b=1:nBoot
            samp = xvec(randsample(n,n,true));
            boot_d(b) = mean(samp,'omitnan') / std(samp,0,'omitnan');
        end

    case 'paired'
        if size(x,2) ~= 2
            error('For paired design, provide Nx2 matrix [pre post] or [before, after].');
        end
        diffs = x(:,2) - x(:,1);
        diffs = diffs(:);
        n = numel(diffs);
        m = mean(diffs,'omitnan');
        s = std(diffs,0,'omitnan');
        d = m / s;
        df = n-1;
        SE = sqrt( (1/n) + (d.^2 ./ (2*n)) );
        boot_d = nan(nBoot,1);
        for b=1:nBoot
            idx = randsample(n,n,true);
            sdiff = diffs(idx);
            boot_d(b) = mean(sdiff,'omitnan') / std(sdiff,0,'omitnan');
        end

    case 'independent'
        if ~iscell(x) || numel(x)~=2
            error('For independent design, supply {x1,x2} as a 1x2 cell array.');
        end
        x1 = x{1}(:); x2 = x{2}(:);
        n1 = numel(x1); n2 = numel(x2);
        m1 = mean(x1,'omitnan'); m2 = mean(x2,'omitnan');
        s1 = std(x1,0,'omitnan'); s2 = std(x2,0,'omitnan');
        % pooled sd
        sp = sqrt(((n1-1)*s1.^2 + (n2-1)*s2.^2) / (n1 + n2 - 2));
        d = (m1 - m2) / sp;
        df = n1 + n2 - 2;
        SE = sqrt( (n1 + n2) / (n1*n2) + (d.^2 / (2*(n1 + n2))) );
        boot_d = nan(nBoot,1);
        for b=1:nBoot
            s1_idx = randsample(n1,n1,true);
            s2_idx = randsample(n2,n2,true);
            sm1 = mean(x1(s1_idx),'omitnan'); sm2 = mean(x2(s2_idx),'omitnan');
            ss1 = std(x1(s1_idx),0,'omitnan'); ss2 = std(x2(s2_idx),0,'omitnan');
            spb = sqrt(((n1-1)*ss1.^2 + (n2-1)*ss2.^2) / (n1 + n2 - 2));
            boot_d(b) = (sm1 - sm2) / spb;
        end
end

% Hedges' g correction J (approx)
J = 1 - (3 ./ (4*df - 1));
g = d .* J;

% analytic CI
z = norminv(1 - alpha/2);
CI_analytic = [d - z*SE, d + z*SE];

% bootstrap CI (percentile)
CI_boot = prctile(boot_d, [100*alpha/2, 100*(1-alpha/2)]);

% assemble output
out.d = d;
out.g = g;
out.J = J;
out.SE_d = SE;
out.CI_analytic = CI_analytic;
out.CI_bootstrap = CI_boot;
out.df = df;
switch typ
    case 'one', out.n = n;
    case 'paired', out.n = n;
    case 'independent', out.n = [n1 n2];
end

if verbose
    fprintf('d = %.3f, g = %.3f, SE = %.3f, 95%% CI analytic [%.3f, %.3f], 95%% CI boot [%.3f, %.3f]\n', ...
        out.d, out.g, out.SE_d, out.CI_analytic(1), out.CI_analytic(2), out.CI_bootstrap(1), out.CI_bootstrap(2));
end

end
