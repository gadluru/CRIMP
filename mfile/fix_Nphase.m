function [Data, para] = fix_Nphase(kSpace, kSpace_info, Image, para)
nset = max(kSpace_info.set(:))+1;
sx = para.Recon.sx;
no_comp = para.Recon.no_comp;
nSMS = para.Recon.nSMS;
cardiac_signal = para.Recon.self_gating.cardiac_signal;
resp_signal = para.Recon.self_gating.resp_signal;

%% pick systole and diastole from cardiac signal
sys = local_max(-cardiac_signal);
dia = local_max(cardiac_signal);
dia(dia<sys(2)) = [];
dia(dia>sys(end-1)) = [];

para.Recon.self_gating.sys = sys;
para.Recon.self_gating.dia = dia;

para.Recon.self_gating.cardiac_signal = cardiac_signal;
para.Recon.self_gating.sys = sys;

%% display self-gating result 
figure
plot(cardiac_signal, 'LineWidth', 2)
hold on
plot(sys,cardiac_signal(sys), '.', 'MarkerSize', 20)
plot(dia,cardiac_signal(dia), '.', 'MarkerSize', 20)
legend('cardiac signal', 'self-gated systole', 'self-gated diastole')
xlabel 'Time Frame'
ylabel 'Cardiac Signal [a. u.]'
set(gca, 'FontSize', 16)

%% 
nor_sl = 6;
nor_one_frame = para.Recon.nor_sl;
para.Recon.nor = nor_one_frame;
ray_idx_start_sys = (sys-1)*nor_sl+1;

ray_idx_start_dia = (dia-1)*nor_sl+1;

Ns2d = round(mean((ray_idx_start_dia - ray_idx_start_sys(2:end-2))/nor_one_frame));
Nd2s = round(mean((ray_idx_start_sys(3:end-1) - ray_idx_start_dia)/nor_one_frame));
Nphase = Ns2d + Nd2s;

Ncycle = length(ray_idx_start_dia);

for i=1:nset
    set_temp = kSpace_info.set==i-1;
    kSpace_temp = kSpace(:,set_temp,:);
    theta_temp = kSpace_info.angle_mod(set_temp);
    phase_temp = kSpace_info.phase_mod(set_temp);
    Image_temp = squeeze(Image(:,:,:,i,:));
    
    kSpace_reorder = zeros(sx,nor_one_frame,Nphase,Ncycle,no_comp);
    theta_reorder = zeros(1,nor_one_frame,Nphase,Ncycle);
    phase_reorder = repmat([0:nSMS-1],[1,nor_one_frame/nSMS*2,Nphase,Ncycle]);
    Image_interp = zeros(sx,sx,Nphase,Ncycle,nSMS);
    resp_signal_reorder = zeros(Nphase,Ncycle);
    
    for icycle=1:Ncycle
        sys_start_temp = ray_idx_start_sys(icycle+1);
        dia_start_temp = ray_idx_start_dia(icycle);
        sys_start_next = ray_idx_start_sys(icycle+2);
        
        % put the systole rays
        kSpace_reorder(:,1:nor_one_frame,1,icycle,:) = kSpace_temp(:,sys_start_temp:sys_start_temp+nor_one_frame-1,:);
        theta_reorder(1,1:nor_one_frame,1,icycle) = theta_temp(sys_start_temp:sys_start_temp+nor_one_frame-1);
        phase_reorder(1,1:nor_one_frame,1,icycle) = phase_temp(sys_start_temp:sys_start_temp+nor_one_frame-1);
        
        % put the systole to diastole rays
        nor_in_between = dia_start_temp - sys_start_temp - nor_one_frame + 1;
        n_frame_in_between = Ns2d-1;
        nor_in_between_one_frame = ceil(nor_in_between/n_frame_in_between/nSMS)*nSMS;
        
        iphase = 2;
        for iframe_in_between = 1:n_frame_in_between
            ray_start_temp = sys_start_temp + nor_one_frame + (iframe_in_between - 1)*nor_in_between_one_frame;
            kSpace_reorder(:,1:nor_in_between_one_frame,iphase,icycle,:) = kSpace_temp(:,ray_start_temp:ray_start_temp+nor_in_between_one_frame-1,:);
            theta_reorder(1,1:nor_in_between_one_frame,iphase,icycle) = theta_temp(ray_start_temp:ray_start_temp+nor_in_between_one_frame-1);
            phase_reorder(1,1:nor_in_between_one_frame,iphase,icycle) = phase_temp(ray_start_temp:ray_start_temp+nor_in_between_one_frame-1);
            iphase = iphase + 1;
        end
        
        % interp image
        frame_temp = sys(icycle+1):dia(icycle);
        Image_temp_for_interp = Image_temp(:,:,frame_temp,:);
        Image_temp_for_interp = permute(Image_temp_for_interp,[3,1,2,4]);
        Image_temp_for_interp = Image_temp_for_interp(:,:);
        Image_temp_for_interp = interp1(1:length(frame_temp),Image_temp_for_interp,1:(length(frame_temp)-1)/(Ns2d-1):length(frame_temp));
        Image_temp_for_interp = reshape(Image_temp_for_interp,[Ns2d,sx,sx,nSMS]);
        Image_temp_for_interp = permute(Image_temp_for_interp,[2,3,1,4]);
        Image_interp(:,:,1:Ns2d,icycle,:) = Image_temp_for_interp;
        
        % interp resp_signal
        resp_signal_temp = resp_signal(frame_temp);
        resp_signal_reorder(1:Ns2d,icycle) = interp1(1:length(frame_temp),resp_signal_temp,1:(length(frame_temp)-1)/(Ns2d-1):length(frame_temp));
        
        % put the diastole rays
        kSpace_reorder(:,1:nor_one_frame,iphase,icycle,:) = kSpace_temp(:,dia_start_temp:dia_start_temp+nor_one_frame-1,:);
        theta_reorder(1,1:nor_one_frame,iphase,icycle) = theta_temp(dia_start_temp:dia_start_temp+nor_one_frame-1);
        phase_reorder(1,1:nor_one_frame,iphase,icycle) = phase_temp(dia_start_temp:dia_start_temp+nor_one_frame-1);
        iphase = iphase + 1;
        
        % put the diastole to systole rays
        nor_in_between = sys_start_next - dia_start_temp - nor_one_frame + 1;
        n_frame_in_between = Nd2s-1;
        nor_in_between_one_frame = ceil(nor_in_between/n_frame_in_between/nSMS)*nSMS;
        
        for iframe_in_between = 1:n_frame_in_between
            ray_start_temp = dia_start_temp + nor_one_frame + (iframe_in_between - 1)*nor_in_between_one_frame;
            kSpace_reorder(:,1:nor_in_between_one_frame,iphase,icycle,:) = kSpace_temp(:,ray_start_temp:ray_start_temp+nor_in_between_one_frame-1,:);
            theta_reorder(1,1:nor_in_between_one_frame,iphase,icycle) = theta_temp(ray_start_temp:ray_start_temp+nor_in_between_one_frame-1);
            phase_reorder(1,1:nor_in_between_one_frame,iphase,icycle) = phase_temp(ray_start_temp:ray_start_temp+nor_in_between_one_frame-1);
            iphase = iphase + 1;
        end
        
        % interp image
        frame_temp = dia(icycle)+1:sys(icycle+2)-1;
        Image_temp_for_interp = Image_temp(:,:,frame_temp,:);
        Image_temp_for_interp = permute(Image_temp_for_interp,[3,1,2,4]);
        Image_temp_for_interp = Image_temp_for_interp(:,:);
        Image_temp_for_interp = interp1(1:length(frame_temp),Image_temp_for_interp,1:(length(frame_temp)-1)/(Nd2s-1):length(frame_temp));
        Image_temp_for_interp = reshape(Image_temp_for_interp,[Nd2s,sx,sx,nSMS]);
        Image_temp_for_interp = permute(Image_temp_for_interp,[2,3,1,4]);
        Image_interp(:,:,Ns2d+1:Nphase,icycle,:) = Image_temp_for_interp;
        
        % interp resp_signal
        resp_signal_temp = resp_signal(frame_temp);
        resp_signal_reorder(Ns2d+1:Nphase,icycle) = interp1(1:length(frame_temp),resp_signal_temp,1:(length(frame_temp)-1)/(Nd2s-1):length(frame_temp));
        
    end
    
    % pre-interpolate
    nor_max = size(kSpace_reorder,2);
    kSpace_reorder = reshape(kSpace_reorder,[sx,nor_max,Nphase*Ncycle,no_comp]);
    theta_reorder = reshape(theta_reorder,[1,nor_max,Nphase*Ncycle]);
    phase_reorder(:,nor_max+1:end,:,:) = [];
    phase_reorder = reshape(phase_reorder,[1,nor_max,Nphase*Ncycle]);
    
    [kx,ky] = get_k_coor(sx,theta_reorder,0,sx/2+1);
    
    % RING trajectory correction
    correction = para.trajectory_correction.*permute([cos(theta_reorder);sin(theta_reorder)],[4,1,2,3]);
    correction = squeeze(sum(correction,2));
    kx = kx - correction(1,:,:);
    ky = ky - correction(2,:,:);
    
    for isms=1:nSMS
        G{isms} = para.G{i,isms};
    end
    Data{i}.kSpace = GROG.SMS_GROG(kSpace_reorder,kx,ky,phase_reorder,para,G);
    Data{i}.first_guess = reshape(Image_interp,[sx,sx,Nphase*Ncycle,1,nSMS]);
    [Data{i},para] = get_Data_SMS(Data{i},para);
    
end

para.Recon.self_gating.resp_signal = resp_signal_reorder(:);
para.Recon.bins = repmat(diag(true(1,Nphase)),[1,Ncycle]);

