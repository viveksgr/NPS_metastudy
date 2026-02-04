function [reducedCells, roiIdx, contrast_sorted] = select_and_sort_rois(interventionsCell, contrast_map, varargin)
% SELECT_AND_SORT_ROIS  Select and reorder ROI columns based on contrast map
%
% [reducedCells, roiIdx, contrast_sorted] = select_and_sort_rois(interventionsCell, contrast_map)
%
% Inputs
%   interventionsCell : 1xK (or Kx1) cell array. Each cell is (nSubjects_i x nROI).
%                       nROI must be the same for all cells.
%   contrast_map      : 1 x nROI (or nROI x 1) numeric vector of t-contrasts.
%
% Optional name/value pairs
%   'threshold'           - scalar (default 3). ROIs > threshold go to the first block.
%   'include_equal_upper' - logical (default false). If true, contrast == threshold goes to upper group.
%   'sortWithinGroup'     - logical (default true). If true, sorts within each group by contrast (desc).
%   'sortDirUpper'        - 'descend' or 'ascend' (default 'descend').
%   'sortDirLower'        - 'descend' or 'ascend' (default 'descend').
%
% Outputs
%   reducedCells    : cell array of same size as interventionsCell. Each cell is nSubjects_i x nROI_new.
%   roiIdx          : 1 x nROI_new integer vector of original ROI indices in the new order.
%   contrast_sorted : 1 x nROI_new contrast_map(roiIdx).
%
% Example:
%   [reducedCells, idx, csorted] = select_and_sort_rois(myCells, contrast_map, 'threshold', 3);
%
% Notes:
%   - ROIs with contrast == threshold are assigned to lower group by default.
%   - If one group is empty the function will still return the other group.
%   - Function does NOT change subject ordering in each cell.

% ---------- parse inputs ----------
p = inputParser;
addParameter(p,'threshold',3,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'include_equal_upper',false,@islogical);
addParameter(p,'sortWithinGroup',true,@islogical);
addParameter(p,'sortDirUpper','descend',@(s) any(strcmp(s,{'descend','ascend'})));
addParameter(p,'sortDirLower','descend',@(s) any(strcmp(s,{'descend','ascend'})));
parse(p,varargin{:});
th = p.Results.threshold;
inclEq = p.Results.include_equal_upper;
sortWithin = p.Results.sortWithinGroup;
sortU = p.Results.sortDirUpper;
sortL = p.Results.sortDirLower;

% ---------- basic checks ----------
if ~iscell(interventionsCell)
    error('interventionsCell must be a cell array.');
end

% find nROI from first non-empty cell
nROI = [];
for ii=1:numel(interventionsCell)
    if ~isempty(interventionsCell{ii})
        [~, nROI] = size(interventionsCell{ii});
        break;
    end
end
if isempty(nROI)
    error('All cells appear empty.');
end

contrast_map = contrast_map(:)'; % row
if numel(contrast_map) ~= nROI
    error('contrast_map length (%d) does not match nROI (%d).', numel(contrast_map), nROI);
end

% ---------- determine groups ----------
upperIdx = find(contrast_map > th);
lowerIdx = find(contrast_map < th);
eqIdx = find(contrast_map == th);

if inclEq
    upperIdx = [upperIdx, eqIdx];
else
    lowerIdx = [lowerIdx, eqIdx];
end

% ---------- sort within groups ----------
if sortWithin
    if ~isempty(upperIdx)
        [~, ordU] = sort(contrast_map(upperIdx), sortU);
        upperSorted = upperIdx(ordU);
    else
        upperSorted = [];
    end
    if ~isempty(lowerIdx)
        [~, ordL] = sort(contrast_map(lowerIdx), sortL);
        lowerSorted = lowerIdx(ordL);
    else
        lowerSorted = [];
    end
else
    upperSorted = upperIdx;
    lowerSorted = lowerIdx;
end

% final index order
roiIdx = [upperSorted, lowerSorted];
contrast_sorted = contrast_map(roiIdx);

% ---------- build reduced cell array ----------
reducedCells = cell(size(interventionsCell));
for ii = 1:numel(interventionsCell)
    Ci = interventionsCell{ii};
    if isempty(Ci)
        reducedCells{ii} = [];
        continue;
    end
    % verify number of columns
    if size(Ci,2) ~= nROI
        error('Cell %d has %d columns but expected %d.', ii, size(Ci,2), nROI);
    end
    reducedCells{ii} = Ci(:, roiIdx);
end

end
