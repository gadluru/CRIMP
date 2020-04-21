function para = prepare_para(para)
fprintf('Loading parameters...');tic

matObj = matfile([para.dir.load_kSpace_dir,para.dir.load_kSpace_name]);
varlist = who(matObj,'kSpace_info');
if ~isempty(varlist)
    load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name],'kSpace_info')
    if isfield(kSpace_info,'phase_mod')
        para.phase_mod = kSpace_info.phase_mod;
    end
    if isfield(kSpace_info,'angle_mod')
        para.angle_mod = kSpace_info.angle_mod;
    end
    if isfield(kSpace_info,'ray_mod')
        ray_mod = kSpace_info.ray_mod(:,1);
        ray_mod = reshape(ray_mod,kSpace_info.Nr,kSpace_info.Nf);
        ray_mod = ray_mod(1,:);
        para.Recon.PD_frames = ray_mod == 0;
        if sum(para.Recon.PD_frames) ~= length(ray_mod) && ~isfield(para.Recon,'RF_frames')
            para.Recon.RF_frames = ray_mod == ray_mod(sum(para.Recon.PD_frames)+1);
            para.Recon.All_rays = para.Recon.PD_frames + para.Recon.RF_frames;
        end
    end
    if isfield(kSpace_info,'line_type')
        ray_mod = kSpace_info.line_type;
        para.Recon.PD_frames = false(1,max(kSpace_info.frames));
        para.Recon.PD_frames(1:kSpace_info.frames(find(diff(ray_mod==0)))) = true;
        para.Recon.RF_frames = ~para.Recon.PD_frames;
    end
end

para.time = datestr(clock,'yymmdd_hhMMSS');


switch para.Recon.interp_method
    case 'NUFFT'
    para.ifNUFFT = 1;
end

if para.phase_mod == 0
    nSMS = 1;
else
    nSMS = 3;
end
para.Recon.nSMS = nSMS;

disp('RawData:')
disp([para.dir.load_kSpace_dir para.dir.load_kSpace_name])

para.dir.save_recon_img_name= sprintf('RECON_%s_%s.mat',para.dir.load_kSpace_name(1:end-4), para.time);

if isempty(dir(para.dir.save_recon_img_mat_dir))
    mkdir(para.dir.save_recon_img_mat_dir);
end

if isempty(dir([pwd,'/RawData/']))
    mkdir([pwd,'/RawData/']);
end

kSpace_data_dir  = para.dir.load_kSpace_dir;
kSpace_data_name = para.dir.load_kSpace_name;

PCA_name = strcat(kSpace_data_name(1:end-4),'_PCA.mat');

PCA_dir_oringinal = [kSpace_data_dir,PCA_name];
if isempty(dir(PCA_dir_oringinal))
    para.dir.PCA_dir = [pwd,'/RawData/',PCA_name];
else
    para.dir.PCA_dir = PCA_dir_oringinal;
end

PCA_precomputed = dir(para.dir.PCA_dir);
para.Recon.ifPCA = ~isempty(PCA_precomputed);

matObj = matfile(para.dir.PCA_dir);
varlist = who(matObj,'kSpace_info');
if ~isempty(varlist)
    load(para.dir.PCA_dir,'kSpace_info')
    if isfield(kSpace_info,'phase_mod')
        para.phase_mod = kSpace_info.phase_mod;
    end
    if isfield(kSpace_info,'angle_mod')
        para.angle_mod = kSpace_info.angle_mod;
    end
    if isfield(kSpace_info,'ray_mod')
        ray_mod = kSpace_info.ray_mod(:,1);
        ray_mod = reshape(ray_mod,kSpace_info.Nr,kSpace_info.Nf);
        ray_mod = ray_mod(1,:);
        para.Recon.PD_frames = ray_mod == 0;
        if sum(para.Recon.PD_frames) ~= length(ray_mod)
            para.Recon.RF_frames = ray_mod == ray_mod(sum(para.Recon.PD_frames)+1);
            para.Recon.All_rays = para.Recon.PD_frames + para.Recon.RF_frames;
        end
    end
    if isfield(kSpace_info,'line_type')
        ray_mod = kSpace_info.line_type;
        if isfield(kSpace_info,'frames')
            para.Recon.PD_frames = false(1,max(kSpace_info.frames));
            para.Recon.PD_frames(1:kSpace_info.frames(find(diff(ray_mod==0)))) = true;
        else
            para.Recon.PD_frames = false(1,length(kSpace_info.MeasurementTime_seconds));
            para.Recon.PD_frames(1:find(diff(ray_mod))/kSpace_info.RadialViews) = true;
        end
        para.Recon.RF_frames = ~para.Recon.PD_frames;
    end
end

if ~isfield(para.Recon,'PD_frames') && exist('kSpace_info')
    if isfield(kSpace_info,'MeasurementTime_seconds')
        para.Recon.PD_frames = false(1,length(kSpace_info.MeasurementTime_seconds));
        para.Recon.PD_frames(1:kSpace_info.ProtonDensityScans) = true;
        para.Recon.RF_frames = ~para.Recon.PD_frames;
    end
end

if exist('kSpace_info')
    para.kSpace_info = kSpace_info;
end

if isempty(dir(para.dir.save_recon_img_mat_dir))
    mkdir(para.dir.save_recon_img_mat_dir)
end

para.CPUtime.load_para_time = toc;toc;fprintf('\n');
end
