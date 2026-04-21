function [meanOfRowMeans, stdOfRowMeans, rowMeans,Xsel] = countThresholdedROIs(X, argpos)
% COUNTTHRESHOLDEDROIS Reduce to selected columns, then summarize row means.
%
% USAGE:
%   [meanOfRowMeans, stdOfRowMeans, rowMeans] = countThresholdedROIs(X, argpos)
%
% INPUTS:
%   X       - [nSubjects x nROI] numeric matrix
%   argpos  - logical vector selecting ROI columns to retain
%
% OUTPUTS:
%   meanOfRowMeans - mean across subjects of the row-wise column means
%   stdOfRowMeans  - standard err across subjects of the row-wise column means
%   rowMeans       - [nSubjects x 1] mean across selected columns for each subject
%
% NOTES:
%   - NaNs are ignored when averaging within rows and across rows.

if nargin < 2
    error('Need X and argpos.');
end

if ~ismatrix(X) || ~isnumeric(X)
    error('X must be a numeric 2D matrix.');
end

if ~islogical(argpos)
    error('argpos must be a logical vector.');
end

argpos = argpos(:);
if size(X, 2) ~= numel(argpos)
    error('Length of argpos must match the number of columns in X.');
end

Xsel = X(:, argpos);
rowMeans = mean(Xsel, 2, 'omitnan');

meanOfRowMeans = mean(rowMeans, 'omitnan');
stdOfRowMeans = std(rowMeans, 0, 'omitnan')./sqrt(length(rowMeans));

end
