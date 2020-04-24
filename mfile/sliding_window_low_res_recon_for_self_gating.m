function [Image, para] = sliding_window_low_res_recon_for_self_gating(kSpace,para)
% [Image, para] = sliding_window_low_res_recon_for_self_gating(kSpace, para)
% Perform a preliminary reconstruction for self-gating

%% pre-interpolation
sx = size(kSpace,1);
no_comp = size(kSpace,3);
nset = max(para.kSpace_info.set(:))+1;
nor_sl = 6;

for i=1:nset
    set = para.kSpace_info.set==i-1;

    kSpace_radial = kSpace(:,set,:);
    theta = para.kSpace_info.angle_mod(set);
    phase = para.kSpace_info.phase_mod(set);
    
    nor_one_frame = para.Recon.nor_sl;
    nor_total = size(kSpace_radial,2);
    
    nof = floor((nor_total-(nor_one_frame-nor_sl))/nor_sl);
    nor_total = nof*nor_sl+(nor_one_frame-nor_sl);
    
    kSpace_radial(:,nor_total+1:end,:) = [];
    theta(nor_total+1:end) = [];
    phase(nor_total+1:end) = [];
    
    kSpace_sl = zeros(sx,nor_one_frame,nof,no_comp);
    theta_sl = zeros(1,nor_one_frame,nof);
    phase_sl = zeros(1,nor_one_frame,nof);
    
    for iframe=1:nof
        ray_idx = (iframe-1)*nor_sl+1:(iframe-1)*nor_sl+nor_one_frame;
        kSpace_sl(:,:,iframe,:) = kSpace_radial(:,ray_idx,:);
        theta_sl(1,:,iframe) = theta(ray_idx);
        phase_sl(1,:,iframe) = phase(ray_idx);
    end
    
    [kx,ky] = get_k_coor(sx,theta_sl,0,sx/2+1);

    correction = para.trajectory_correction.*permute([cos(theta_sl);sin(theta_sl)],[4,1,2,3]);
    correction = squeeze(sum(correction,2));
    kx = kx - correction(1,:,:);
    ky = ky - correction(2,:,:);

    for j=1:para.Recon.nSMS
        G{j} = para.G{i,j};
    end
    [Data{i}.kSpace,Data{i}.G] = GROG.SMS_GROG(kSpace_sl,kx,ky,phase_sl,para,G);
    para.Recon.nor = nor_one_frame;
    [Data{i},para] = get_Data_SMS(Data{i},para);

    para.Recon.no_comp = no_comp;
end

para.kSpace_info.TimeStamp = para.kSpace_info.TimeStamp(1:nset*nor_sl:end);

siz = size(Data{1}.first_est);
siz(4) = nset;
Image = zeros(siz,'single');

%% recon reference image
para.Recon.type = 'seperate SMS test';
para.Recon.noi = 50;

for i=1:nset
    scale = max(abs(Data{i}.first_est(:)));
    para.Recon.weight_tTV = scale*para.weight_tTV;
    para.Recon.weight_sTV = scale*para.weight_sTV;
    for j=1:para.Recon.no_comp
        k_temp = Data{i}.kSpace(:,:,:,j,:,:,:);
        kSpace_temp(:,j) = k_temp(Data{i}.mask);
    end
    Data{i}.kSpace = kSpace_temp;
    Image(:,:,:,i,:) = STCR_conjugate_gradient(Data{i},para);
    clear kSpace_temp
end
