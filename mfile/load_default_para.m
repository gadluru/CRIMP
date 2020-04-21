function para = load_default_para(option)
if ~exist('option')
    option = '';
end
if contains(option,'GROG')
    para.Recon.interp_method = 'GROG';
    para.over_sampling = 1;
    para.core_size = [1,1];
end
if contains(option,'NUFFT')
    para.Recon.interp_method = 'NUFFT';
    para.over_sampling = 1.5;
    para.core_size = [6,6];
end
if contains(option,'GPU')
    para.setting.ifGPU = 1;
else
    para.setting.ifGPU = 0;
end

para.setting.save_frequency = 20000;
para.setting.debug = 1;
para.setting.plot = 1;

para.beta_square = eps('single');

para.weight_sTV = 0.001;
para.weight_tTV = 0.02;

para.step_size = 2;
para.Recon.crop_half_FOV = 1;
para.image_orintation = 7;
para.noi = 150;
para.Recon.break = 0;

para.dir.save_recon_img_mat_dir = strcat(pwd,'/ReconData/mat_files/');