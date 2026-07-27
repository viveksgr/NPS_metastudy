% MAKE_TEST2_FIGS  Test-2 decoding bar plots.
% Bar 1: overall balanced accuracy from full z, readout axis alone, null space alone.
% Bar 2: same, broken out per intervention (per-group recall).
% Error bars = bootstrap over studies; stars = vs permutation floor.
setup;
S = load('datamat.mat'); cl_mat = S.cl_mat;

% COORDS selectable: 'c3' -> {C1,C2,C3}, 'c5' -> add NPS,SIIPS.
% Set at top or pass make_test2_figs from a wrapper; defaults to c3.
if ~exist('WHICH','var'), WHICH = 'c3'; end
switch lower(WHICH)
    case 'c5'
        COORDS = {'C1','C2','C3','NPS','SIIPS'}; tag = 'c5'; ttlSuf = ' (5D: C1-3+NPS+SIIPS)';
    otherwise
        COORDS = {'C1','C2','C3'};               tag = 'c3'; ttlSuf = ' (3D: C1-3)';
end
matf = sprintf('outputs/test2_R_%s.mat', tag);

% cache so the figures can be restyled without re-running the permutations
if exist(matf,'file')
    load(matf,'R2');
else
    R2 = run_test2('COORDS',COORDS,'N_PERM',500,'N_BOOT',1000);
    save(matf,'R2');
end

plotTest2Acc(R2, cl_mat, ['Test 2 - decoding by subspace (LOSO)' ttlSuf]);
exportgraphics(gcf, sprintf('outputs/fig_test2_acc_overall_%s.png',tag), 'Resolution', 200);

plotTest2Acc(R2, cl_mat, ['Test 2 - decoding by subspace, per intervention' ttlSuf], 'byGroup', true);
exportgraphics(gcf, sprintf('outputs/fig_test2_acc_bygroup_%s.png',tag), 'Resolution', 200);

writetable(array2table([R2.accG R2.pG], 'RowNames', R2.GRP, ...
    'VariableNames', {'full','readout','null','p_full','p_readout','p_null'}), ...
    sprintf('outputs/test2_bygroup_%s.csv',tag), 'WriteRowNames', true);
fprintf('test2 %s figures + csv written\n', tag);
