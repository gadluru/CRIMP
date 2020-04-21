addpath ./mfile/
ccc
%%
% This is the script runs the reconstruction described in the CRIMP paper
% The reconstruction deals with images of the size >2Gb, so the processing
% is slow. The entire processing requires ~20Gb memory. The processing on a
% GPU is >3 hours, and can be slower when running on CPU. The
% reconstruction can be partially tested by settings in the function
% 'recon_2D_ungated_SPGR_mSMS_low_rank_car_ph_resolved_ADMM_auto' line
% 77-85.

cd ./reconstruction/
%% sort k-space data together
% GitHub allows a single file <100mb. The k-space was seperated into 3
% files. This script will combine those into a single *.mat file.
run('./combine_kspace.m')

%% the main reconstruction script
run('./recon_CRIMP_ADMM.m')
% reconstruction results will be saved at:
% ./reconstruction/ReconData/mat_files/
