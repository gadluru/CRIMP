function [SI_one_slice,Mz] = simu_steady_state_with_motion(flip_angle_slice,Mz,M0,TR,Slab_T1,max_disp,heart_rate)
T = 1:3000; % number of alpha pulses
total_time = T(end)*TR/1000;
Ncycle = total_time/60*heart_rate;
max_disp = 100*max_disp; % one slice thickness = 200 pixels
Motion = round(sin(T/T(end)*Ncycle*2*pi)*max_disp);

interleaving_slices = repmat([1,2,3],[1,length(T)/3]);
SI_one_slice = zeros(9, 3000, 'single');

for i=T
    flip_angle_all = zeros(3800, 1, 'single');
    
    SMS_slices = [201:1200,801:1800,1401:2400]+400;
    if interleaving_slices(i) == 2
        SMS_slices = SMS_slices + 200;
    elseif interleaving_slices(i) == 3
        SMS_slices = SMS_slices + 400;
    end
    SMS_slices = SMS_slices+Motion(i);
    
    flip_angle_all(SMS_slices(1:1000)) = flip_angle_all(SMS_slices(1:1000)) + flip_angle_slice;
    flip_angle_all(SMS_slices(1001:2000)) = flip_angle_all(SMS_slices(1001:2000)) + flip_angle_slice;
    flip_angle_all(SMS_slices(2001:3000)) = flip_angle_all(SMS_slices(2001:3000)) + flip_angle_slice;
    
    Mz(:,end+1) = Mz(:,end).*cos(flip_angle_all/180*pi);
%     SI(:,i) = Mz(:,end-1).*sin(flip_angle_all/180*pi); % this is a
%     constant so not so important
    for iii=1:9
        slice_location = (201:1200) + 200*(iii-1) + Motion(i) + 400;
        SI_one_slice(iii,i) = sum(Mz(slice_location,end-1).*sin(flip_angle_slice/180*pi));
    end
    Mz(:,end+1) = M0 - (M0 - Mz(:,end)).*exp(-TR./Slab_T1);
end

SI_one_slice_temp([1,4,7],:) = SI_one_slice([1,4,7],1:3:end);
SI_one_slice_temp([2,5,8],:) = SI_one_slice([2,5,8],2:3:end);
SI_one_slice_temp([3,6,9],:) = SI_one_slice([3,6,9],3:3:end);

SI_one_slice = SI_one_slice_temp;
end