addpath ./mfile/
ccc
%%
% This is the script to run CNR/SNR comparision between saturation recovery
% sequence and SPGR sequences (CRIMP), and reproduce Figure 5, 6 in paper.

cd ./SNR_CNR_simulation/

%% repoduce Figure 5a
run('Simulation_SR_SPGR_slice_group_cnr.m')
fprintf('press key to continue...\n')
pause
%% repoduce Figure 6a
run('Simulation_SR_SPGR_FA_cnr.m')
fprintf('press key to continue...\n')
pause
%% repoduce Figure 6a
run('Simulation_SR_SPGR_slice_group_snr.m')
fprintf('press key to continue...\n')
pause
%% repoduce Figure 6b
run('Simulation_SR_SPGR_FA_snr.m')
fprintf('press key to continue...\n')
pause

cd ../