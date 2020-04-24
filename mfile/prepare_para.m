function para = prepare_para(para)
fprintf('Loading parameters...');tic

matObj = matfile([para.dir.load_kSpace_dir,para.dir.load_kSpace_name]);
varlist = who(matObj,'kSpace_info');
if ~isempty(varlist)
    load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name],'kSpace_info')
    para.Recon.nSMS = max(kSpace_info.phase_mod) + 1;
    para.kSpace_info = kSpace_info;
end

para.time = datestr(clock,'yymmdd_hhMMSS');

disp('RawData:')
disp([para.dir.load_kSpace_dir, para.dir.load_kSpace_name])

para.dir.save_recon_img_name= sprintf('RECON_%s_%s.mat',para.dir.load_kSpace_name(1:end-4), para.time);

if isempty(dir(para.dir.save_recon_img_dir))
    mkdir(para.dir.save_recon_img_dir);
end

if isempty(dir([pwd,'/RawData/']))
    mkdir([pwd,'/RawData/']);
end

para.CPUtime.load_para_time = toc;toc;fprintf('\n');
end
