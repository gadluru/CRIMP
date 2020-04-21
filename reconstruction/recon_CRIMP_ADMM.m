addpath ../mfile/
ccc

all_mat = dir('./RawData/*Ungated2D*PCA*.mat');
file = all_mat(1);

%% set reconstruction parameters
para = load_default_para('GROG');

para.dir.load_kSpace_name = file.name;
para.dir.load_kSpace_dir = [file.folder,'/'];
para.Recon.noi = 30; % number of iterations. This is changed several times during the reconstruction.
para.kSpace_center = 145; % k-space center point along readout
para.weight_sTV = 0.00; % spatial constraint weight
para.weight_tTV = 0.04; % temporal constraint weight
para.nor_sl = 12; % number of sliding window rays
para.load_flag = 0;
para.setting.ifGPU = 1; % you need a Nvidia GPU memory of ~24Gb, otherwise set to 0
para.setting.plot = 1; % show process during reconstruction

para = prepare_para(para);

%% run the main reconstruction function
recon_2D_ungated_SPGR_mSMS_low_rank_car_ph_resolved_ADMM_auto(para);
