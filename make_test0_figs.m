% MAKE_TEST0_FIGS  Per-group Test-0 bar plots: random CV vs LOSO.
setup;
S = load('datamat.mat');
cl_mat = S.cl_mat;

Rcv   = run_test0_bygroup('MODE','cv');
Rloso = run_test0_bygroup('MODE','loso');

plotTest0ByGroup(Rcv,   cl_mat, 'Test 0 by group - random CV (ceiling pass)');
exportgraphics(gcf, 'outputs/fig_test0_bygroup_cv.png', 'Resolution', 200);

plotTest0ByGroup(Rloso, cl_mat, 'Test 0 by group - leave-one-study-out');
exportgraphics(gcf, 'outputs/fig_test0_bygroup_loso.png', 'Resolution', 200);

writetable(Rcv,   'outputs/test0_bygroup_cv.csv');
writetable(Rloso, 'outputs/test0_bygroup_loso.csv');
disp('figures + csv written');
