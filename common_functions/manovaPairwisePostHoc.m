function postStats = manovaPairwisePostHoc(Y, group, varargin)
% manovaPairwisePostHoc  Pairwise Hotelling's T^2 post-hoc for manova1.
%
%   postStats = manovaPairwisePostHoc(Y, group)
%   postStats = manovaPairwisePostHoc(Y, group, 'NBoot',2000, 'Adjust','holm')
%
% INPUTS
%   Y      - n x p matrix of dependent variables (p = 2 in your case)
%   group  - n x 1 grouping variable (numeric, cell, string, categorical)
%
% NAME/VALUE
%   'NBoot'  - bootstrap reps for SE/CI on Mahalanobis distance (default 2000)
%   'Adjust' - p-value adjustment: 'none' (default), 'bonferroni', 'holm'
%   'Alpha'  - CI level (default 0.05 -> 95% CI)
%
% OUTPUT (compatible with plotFixedEffects)
%   postStats.Name     - cell array of pair labels "A vs B"
%   postStats.Estimate - Mahalanobis distance between centroids
%   postStats.SE       - bootstrap SE of distance
%   postStats.Lower    - lower CI bound
%   postStats.Upper    - upper CI bound
%   postStats.pValue   - Hotelling's T^2 p-value (F-approximation), adjusted

    ip = inputParser;
    addParameter(ip,'NBoot',2000,@(x)isnumeric(x)&&isscalar(x));
    addParameter(ip,'Adjust','none', ...
        @(s) any(strcmpi(s,{'none','bonferroni','holm'})));
    addParameter(ip,'Alpha',0.05,@(x)isnumeric(x)&&isscalar(x));
    parse(ip,varargin{:});
    nBoot = ip.Results.NBoot;
    adj   = lower(ip.Results.Adjust);
    alpha = ip.Results.Alpha;

    if istable(group), group = group{:,1}; end
    if isnumeric(group)
        g = categorical(cellstr(num2str(group(:))));
    else
        g = categorical(cellstr(string(group(:))));
    end
    levels = categories(g);
    G = numel(levels);
    p = size(Y,2);

    pairs = nchoosek(1:G, 2);
    nP = size(pairs,1);

    Name     = cell(nP,1);
    Estimate = zeros(nP,1);
    SE       = zeros(nP,1);
    Lower    = zeros(nP,1);
    Upper    = zeros(nP,1);
    pValue   = zeros(nP,1);

    for k = 1:nP
        i = pairs(k,1); j = pairs(k,2);
        Yi = Y(g == levels{i}, :);
        Yj = Y(g == levels{j}, :);
        Yi = Yi(all(isfinite(Yi),2), :);
        Yj = Yj(all(isfinite(Yj),2), :);
        ni = size(Yi,1); nj = size(Yj,1);

        [D, pVal] = hotellingPair(Yi, Yj);

        % bootstrap SE / percentile CI on D
        Dboot = nan(nBoot,1);
        for b = 1:nBoot
            bi = Yi(randi(ni, ni, 1), :);
            bj = Yj(randi(nj, nj, 1), :);
            Dboot(b) = mahalDist(bi, bj);
        end
        Dboot = Dboot(isfinite(Dboot));

        Name{k}     = sprintf('%s vs %s', levels{i}, levels{j});
        Estimate(k) = D;
        SE(k)       = std(Dboot);
        Lower(k)    = quantile(Dboot, alpha/2);
        Upper(k)    = quantile(Dboot, 1 - alpha/2);
        pValue(k)   = pVal;
    end

    % multiple-comparison adjustment
    switch adj
        case 'bonferroni'
            pValue = min(pValue * nP, 1);
        case 'holm'
            [psort, ord] = sort(pValue);
            adjP = min(1, psort .* (nP - (0:nP-1)'));
            adjP = cummax(adjP);
            pValue(ord) = adjP;
    end

    postStats.Name     = Name;
    postStats.Estimate = Estimate;
    postStats.SE       = SE;
    postStats.Lower    = Lower;
    postStats.Upper    = Upper;
    postStats.pValue   = pValue;
end

function D = mahalDist(A, B)
% pooled-covariance Mahalanobis distance between centroids of A and B
    nA = size(A,1); nB = size(B,1);
    mA = mean(A,1); mB = mean(B,1);
    Sp = ((nA-1)*cov(A) + (nB-1)*cov(B)) / (nA + nB - 2);
    diff = (mA - mB)';
    D = sqrt(diff' / Sp * diff);
end

function [D, pVal] = hotellingPair(A, B)
% two-sample Hotelling's T^2 with F approximation
    nA = size(A,1); nB = size(B,1); p = size(A,2);
    mA = mean(A,1); mB = mean(B,1);
    Sp = ((nA-1)*cov(A) + (nB-1)*cov(B)) / (nA + nB - 2);
    diff = (mA - mB)';
    T2 = (nA*nB/(nA+nB)) * (diff' / Sp * diff);
    df1 = p;
    df2 = nA + nB - p - 1;
    F   = T2 * df2 / (p * (nA + nB - 2));
    pVal = 1 - fcdf(F, df1, df2);
    D    = sqrt(diff' / Sp * diff);   % Mahalanobis distance between centroids
end
