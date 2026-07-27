% MAKE_TEST1_FIGS  Test-1 separation-split bar plots.
% Bar 1: mean on-axis (a) vs mean null-space (b) squared separation, SE across pairs.
% Bar 2: same, averaged within each intervention (all pairs containing it).
% Dashed red line = isotropy reference b_iso = (k-1)*a.
setup;
S = load('datamat.mat'); cl_mat = S.cl_mat;

R = run_test1('MODE','all','COORDS',{'C1','C2','C3'},'N_PERM',2000,'N_BOOT',1000);

plotTest1AB(R, cl_mat, 'Test 1 - on-axis vs null-space separation (all pairs)');
exportgraphics(gcf, 'outputs/fig_test1_ab_overall.png', 'Resolution', 200);

plotTest1AB(R, cl_mat, 'Test 1 - separation split by intervention', 'byGroup', true);
exportgraphics(gcf, 'outputs/fig_test1_ab_bygroup.png', 'Resolution', 200);

disp('test1 figures written');
