function BigData_norm = l2normalize_columns(BigData,method)
% L2NORMALIZE_COLUMNS  L2-normalize each column, excluding NaNs
%
%   BigData_norm = l2normalize_columns(BigData)
%
% INPUT:
%   BigData      – [N×M] data matrix (may contain NaNs)
%
% OUTPUT:
%   BigData_norm – [N×M] same size, each column scaled to unit L2 norm
%                   (NaNs are left in place)

[N, M] = size(BigData);
BigData_norm = nan(N, M);

for j = 1:M
    col    = BigData(:, j);
    valid  = ~isnan(col);
    if any(valid)
        % compute L2 norm over non-NaN entries
        norm_j = sqrt( sum(col(valid).^2) );
        if norm_j > 0
            % divide only the valid entries
            if strcmp(method,'L2')
                BigData_norm(valid, j) = col(valid) / norm_j;
            elseif strcmp(method,'InvVar_n')
                BigData_norm(valid, j) = col(valid) / (var(col(valid))/sum(valid));
            elseif strcmp(method,'InvVar')
                BigData_norm(valid, j) = col(valid) / (var(col(valid)));
            else 
                error('Wrong method')
            end

        else
            % if all zeros, leave as zeros (or as-is)
            BigData_norm(valid, j) = 0;
        end
    end
    % BigData_norm(~valid,j) remains NaN
end
end
