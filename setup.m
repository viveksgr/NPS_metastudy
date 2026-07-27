% SETUP  Put the analysis functions and datasets on the MATLAB path.
% Run once per session before calling run_* / make_* from the repo root.
%   >> setup
% load('datamat.mat') then resolves via the path; figures/caches go to outputs/.
addpath('common_functions');   % run_* analyses + plot* helpers
addpath('data');               % datamat.mat, datamat_old.mat (load searches path)
if ~exist('outputs','dir'), mkdir('outputs'); end
