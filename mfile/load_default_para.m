function para = load_default_para(option)
if ~exist('option')
    option = '';
end

if contains(option,'GPU')
    para.setting.ifGPU = 1;
else
    para.setting.ifGPU = 0;
end
para.setting.ifplot = 1;

para.weight_sTV = 0.001;
para.weight_tTV = 0.02;

para.Recon.epsilon = eps('single');
para.Recon.step_size = 2;
para.Recon.noi = 150;
para.Recon.break = 0;

para.dir.save_recon_img_dir = strcat(pwd,'/ReconData/mat_files/');