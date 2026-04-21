function cellmat = conv2cell(mat)

cellmat = cell(1, size(mat,2));
for zz = 1:size(mat,2)
    thr = isnan(mat(:,zz));
    thr_in = diff((thr));
    idx = find(thr_in==1);
    if isempty(idx); idx= size(mat,1); end
    cellmat{zz} = mat(1:idx(end),zz);
end