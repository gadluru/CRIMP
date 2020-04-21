function [Image_PD,Data] = recon_PD_multiple_frames(kSpace_all,kSpace_info,para)

kSpace_info.angle_mod = kSpace_info.angle_mod(:,1:kSpace_info.PD_rays,:);
kSpace_info.phase_mod = kSpace_info.phase_mod(:,1:kSpace_info.PD_rays,:);
kSpace_info.set = kSpace_info.set(:,1:kSpace_info.PD_rays,:);

para.kSpace_info = kSpace_info;

[sx,nor_all,no_comp] = size(kSpace_all);

para.Recon.nSMS = max(kSpace_info.phase_mod(:))+1;
kCenter = para.kSpace_center;
para.setting.plot = 1;

nset = max(kSpace_info.set(:))+1;

%% recon reference 
for i=1:nset
    set = kSpace_info.set==i-1;
    kSpace_radial = kSpace_all(:,set,:);
    theta = kSpace_info.angle_mod(set);
    phase = kSpace_info.phase_mod(set);

    nor_sl = para.nor_sl;
    
    nor_total = size(kSpace_radial,2);
    nof = floor(nor_total/nor_sl);
    nor_total = nof*nor_sl;
    kSpace_radial(:,nor_total+1:end,:) = [];
    theta(nor_total+1:end) = [];
    phase(nor_total+1:end) = [];
    
    kSpace_radial = reshape(kSpace_radial,[sx,nor_sl,nof,no_comp]);
    theta = reshape(theta,[1,nor_sl,nof]);
    phase = reshape(phase,[1,nor_sl,nof]);
    
    [kx,ky] = get_k_coor(sx,theta,0,kCenter);

    correction = para.trajectory_correction.*permute([cos(theta);sin(theta)],[4,1,2,3]);
    correction = squeeze(sum(correction,2));
    kx = kx - correction(1,:,:);
    ky = ky - correction(2,:,:);
    %phase = get_phase(phase);
    
    [Data{i}.kSpace,Data{i}.G] = GROG.SMS_GROG(kSpace_radial,kx,ky,phase,para);
    para.Recon.nor = nor_sl;
    [Data{i},para] = get_Data_SMS(Data{i},para);

    para.Recon.no_comp = no_comp;
end

siz = size(Data{1}.first_est);
siz(4) = nset;
Image_PD = zeros(siz,'single');
% para.setting.ifGPU = 0;

for i=1:nset
    scale = max(abs(Data{i}.first_est(:)));
    para.Recon.weight_tTV = scale*0.05;
    para.Recon.weight_sTV = scale*0.008;
    para.setting.plot = 0;
    para.Recon.type = 'seperate SMS test';
    para.Recon.noi = 150;
    for j=1:para.Recon.no_comp
        k_temp = Data{i}.kSpace(:,:,:,j,:,:,:);
        kSpace(:,j) = k_temp(Data{i}.mask);
    end
    Data{i}.kSpace = kSpace;
    Image_PD(:,:,:,i,:) = STCR_conjugate_gradient_low_rank_bins(Data{i},para);
    clear kSpace
end

Image_PD = abs(crop_half_FOV(Image_PD(:,:,:,:)));
Image_PD = squeeze(Image_PD);
order = vec((1:nset)'+([1,para.Recon.nSMS:-1:2]-1)*nset);
Image_PD = Image_PD(:,:,:,order);
