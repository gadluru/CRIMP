addpath ./mfile/
ccc
%%
% This is the script to run motion simulation for Figure 8 and Supporting
% Information Figure S8.

cd ./motion_simulation/

%% repoduce Figure 8
run('plot_figure_8.m')
pause
fprintf('press key to continue...\n')
%% reproduce Figure S8
% The following two scriptes are commonted since long computation time. The
% data is provided in './motion_simulation/data.mat'. Interesting readears
% are encouraged to run those scripts.

% run('simu_motion.m')
% run('make_data_figure_S8.m')
run('plot_figure_S8.m')

cd ../