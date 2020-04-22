function recon_2D_ungated_SPGR_mSMS_low_rank_car_ph_resolved_ADMM_auto(para)
%   recon_2D_ungated_SPGR_mSMS_low_rank_car_ph_resolved_ADMM_auto(kSpace_all,kSpace_info,para)
%
%   This is the main reconstruction function
%
%   Inputs:
%           'para': reconstruction parameters
%               para.dir: directories
%               para.setting: reconstruciton settings
%               para.Recon: useful parameters during reconstruction
%
%   Load data upder 'para.dir':
%           'kSpace': all the acquired radial k-space data in the
%           temporal order. The dimentions are [sx,nor,ncoil], with sx the
%           number of points along each ray, nor the total number of rays
%           and ncoil the number of coils.
%
%           'kSpace_info': include some information for k-space data
%               kSpace_info.angle_mod: radial projection angle
%               kSpace_info.phase_mod: CAIPIRINHA phase modulation pattern
%               kSpace_info.set: interleaving slice groups index
%               kSpace_info.TimeStamp: time stamps
% 
% Copyright: University of Utah Cardiovascular MRI Group
% https://medicine.utah.edu/radiology/radiology-research/research-labs/dibella2/
%
% Contact:
% Ye Tian, phye1988@gmail.com

%% load data and save parameters
load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name], 'kSpace', 'kSpace_info')
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'para','-v7.3')

%% set some parameters
set = kSpace_info.set==0;
set = find(set); set = set(1:600);
nset = max(kSpace_info.set(:))+1;
nSMS = max(kSpace_info.phase_mod(:))+1; % number of SMS slices
para.Recon.nSMS = nSMS;
[sx,nor_all,no_comp] = size(kSpace);

para.Recon.sx = sx;
para.Recon.no_comp = no_comp;

non_steady_state_rays = 360;

%% RING trajectory correction using PD rays
para.trajectory_correction = RING_SMS(kSpace(:,set,:), kSpace_info.angle_mod(set), nSMS);

%% reconstruct proton density (PD) weighted images
Image_PD = recon_PD_multiple_frames(kSpace(:,1:kSpace_info.PD_rays,:),kSpace_info,para);
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image_PD','-append')

%% cut PD rays and non_steaty_state rays
kSpace(:,1:kSpace_info.PD_rays+non_steady_state_rays,:) = [];
kSpace_info.angle_mod(:,1:kSpace_info.PD_rays+non_steady_state_rays,:) = [];
kSpace_info.phase_mod(:,1:kSpace_info.PD_rays+non_steady_state_rays,:) = [];
kSpace_info.set(:,1:kSpace_info.PD_rays+non_steady_state_rays,:) = [];
kSpace_info.TimeStamp(:,1:kSpace_info.PD_rays+non_steady_state_rays,:) = [];
para.kSpace_info = kSpace_info;

%% estimate SMS-GROG operator
fprintf('Estimate SMS-GROG operator...\n')
for iset=1:nset
    tic
    fprintf(sprintf('SMS slice group %g\n', iset))
    set_idx = kSpace_info.set==iset-1;
    kSpace_temp = kSpace(:,set_idx,:);
    theta_temp = kSpace_info.angle_mod(set_idx);
    phase_temp = kSpace_info.phase_mod(set_idx);
    for isms = 1:nSMS
        phase_temp_SMS = phase_temp == nSMS-1;
        kSpace_temp_SMS = kSpace_temp(:,phase_temp_SMS,:);
        kSpace_temp_SMS = permute(kSpace_temp_SMS,[1,2,4,3]);
        theta_temp_SMS = theta_temp(phase_temp_SMS);
        [kx_temp, ky_temp] = get_k_coor(sx,theta_temp_SMS,0,sx/2+1);
        [G{iset,isms}.Gx, G{iset,isms}.Gy] = GROG.get_Gx_Gy(kSpace_temp_SMS, kx_temp, ky_temp);
    end
    toc
end
clear *temp* i*
para.G = G;

% truncate k-space for testing
% kSpace = kSpace(:,1:2016,:);
% para.kSpace_info.set = para.kSpace_info.set(1:2016);
% para.kSpace_info.angle_mod = para.kSpace_info.angle_mod(1:2016);
% para.kSpace_info.phase_mod = para.kSpace_info.phase_mod(1:2016);
% kSpace_info.set = para.kSpace_info.set(1:2016);
% kSpace_info.angle_mod = para.kSpace_info.angle_mod(1:2016);
% kSpace_info.phase_mod = para.kSpace_info.phase_mod(1:2016);
%
%% recon reference image for self-gating
Image = sliding_window_low_res_recon_for_self_gating(kSpace,para);
% save('temp1.mat','Image')
%% slef-gating 
% load('temp1.mat','Image')
[cardiac_signal, resp_signal, para.mask] = get_cardiac_signal_SPGR(Image);
para.Recon.self_gating.resp_signal = resp_signal;
para.Recon.self_gating.cardiac_signal = cardiac_signal;

%% fix the number of cardiac phases for each cardiac cycle
fprintf('Fix number of cardiac phases in all cardiac cycles...\n')
tic
[Data, para] = fix_Nphase(kSpace, kSpace_info, Image, para);
clear Image kSpace
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'para','-append')
toc
%% 2nd reconstruction stage initialization
% put all sets into one Data structure
para.Recon.noi = 50;
Data_all = Data{1};
scale = max(abs(Data{1}.first_est(:)));
for i=2:nset
    Data_all.kSpace = cat(6,Data_all.kSpace,Data{i}.kSpace);
    Data_all.mask = cat(6,Data_all.mask,Data{i}.mask);
    scale = max([scale;abs(Data{i}.first_est(:))]);
    Data_all.sens_map = cat(6,Data_all.sens_map,Data{i}.sens_map);
    Data_all.first_guess = cat(6,Data_all.first_guess,Data{i}.first_guess);
end
clear Data

% sort Cartesian k-space that only have data. This can save huge memory
for i=1:para.Recon.no_comp
    k_temp = Data_all.kSpace(:,:,:,i,:,:,:);
    kSpace_aline(:,i) = k_temp(Data_all.mask);
end
Data_all.kSpace = kSpace_aline;
clear kSpace_aline k_temp

% some reconstruction setting
para.Recon.weight_tTV = scale*0.04;
para.Recon.weight_sTV = scale*0.00;
para.Recon.type = 'seperate SMS test';
para.Recon.noi = 50;
clearvars -except Data_all para
Data_all = rmfield(Data_all,'first_est');

% reorder image slices 
[sx,sy,nof,~,nSMS,nset] = size(Data_all.first_guess);
nslice = nSMS*nset;
order = 1:nSMS:nslice;
for i=nSMS:-1:2
    order = [order,i:nSMS:nslice];
end
[~,order_back] = sort(order);

Data_all.first_guess = reshape(Data_all.first_guess,[sx,sy,nof,nslice]);
Data_all.first_guess = Data_all.first_guess(:,:,:,order);
Data_all.first_guess = single(Data_all.first_guess);

%% estimate rigid translations from diastolic images
fprintf('Start rigid translation estimation...\n')
% figure
% plot(para.Recon.self_gating.resp_signal)
dia_loc = round(size(para.Recon.bins,1)/2+1);
Image_dia = abs(crop_half_FOV(Data_all.first_guess(:,:,para.Recon.bins(dia_loc,:),:)));
Image_dia = permute(Image_dia,[1,2,4,3]);

para.mask = imresize(para.mask,[sx/2,sx/2]);
para.mask = para.mask(:,:,order);

[~,shifts_dia,~] = rigid_reg_3D(Image_dia,para.mask);

% interpolate translations to all cardiac phases
shifts_interp = interp_rigid_shifts(shifts_dia,para.Recon.self_gating.resp_signal,para.Recon.bins(dia_loc,:));
shifts_interp = shifts_interp - mean(shifts_interp,2);

Data_all.shifts = shifts_interp;
para.shifts = Data_all.shifts;

Data_all.size.Nphase = size(para.Recon.bins,1);
Data_all.size.Ncycle = nof/Data_all.size.Nphase;
Data_all.size.order = order;
Data_all.size.order_back = order_back;
Data_all.size.nSMS = nSMS;
Data_all.size.nset = nset;
clearvars -except Data_all para

%% patch tracking
fprintf('Start patch tracking...\n')
patch_size = [5,5,2];
search_size = [9,9,4];
patch_shift = [5,5,2];
Nphase = size(para.Recon.bins,1);
Data_all.llr = Patch_tracking_3D_with_guide_random_search(Data_all.first_guess,Nphase,Data_all.shifts,patch_size,search_size,patch_shift);

para.setting.ifGPU = 0;
para.Recon.noi = 30;
clearvars -except Data_all para
%% ADMM solver
im_bin = LLR_patch_tracking_4D_ADMM(Data_all,para);

%% save reconstructed images
im_bin = permute(im_bin,[1,2,4,3]);
im_bin = abs(crop_half_FOV(im_bin));
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'im_bin','para','-append')

Image_sys = im_bin(:,:,:,para.Recon.bins(1,:));
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image_sys','-append')

dia_loc = round((size(para.Recon.bins,1)+1)/2);
Image_dia = im_bin(:,:,:,para.Recon.bins(dia_loc,:));
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image_dia','-append')

return






      