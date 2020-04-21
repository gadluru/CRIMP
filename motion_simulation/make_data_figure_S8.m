ccc

all_mat = dir('HR*');
N = length(all_mat);
%% calculate min max

for i=1:N
    load(fullfile(all_mat(i).folder,all_mat(i).name),'SI_with_motion')
    if strfind(all_mat(i).name,'60')
        HR = 60;
    elseif strfind(all_mat(i).name,'70')
        HR = 70;
    elseif strfind(all_mat(i).name,'80')
        HR = 80;
    elseif strfind(all_mat(i).name,'90')
        HR = 90;
    elseif strfind(all_mat(i).name,'100')
        HR = 100;
    elseif strfind(all_mat(i).name,'110')
        HR = 110;
    elseif strfind(all_mat(i).name,'120')
        HR = 120;
    end
    if strfind(all_mat(i).name,'_0.5_')
        disp = 0.5;
    elseif strfind(all_mat(i).name,'_0.6_')
        disp = 0.6;
    elseif strfind(all_mat(i).name,'_0.7_')
        disp = 0.7;
    elseif strfind(all_mat(i).name,'_0.8_')
        disp = 0.8;
    elseif strfind(all_mat(i).name,'_0.9_')
        disp = 0.9;
    elseif strfind(all_mat(i).name,'_1_')
        disp = 1;
    elseif strfind(all_mat(i).name,'_1.1_')
        disp = 1.1;
    elseif strfind(all_mat(i).name,'_1.2_')
        disp = 1.2;
    elseif strfind(all_mat(i).name,'_1.3_')
        disp = 1.3;
    elseif strfind(all_mat(i).name,'_1.4_')
        disp = 1.4;
    elseif strfind(all_mat(i).name,'_1.5_')
        disp = 1.5;
    elseif strfind(all_mat(i).name,'_1.6_')
        disp = 1.6;
    elseif strfind(all_mat(i).name,'_1.7_')
        disp = 1.7;
    elseif strfind(all_mat(i).name,'_1.8_')
        disp = 1.8;
    elseif strfind(all_mat(i).name,'_1.9_')
        disp = 1.9;
    elseif strfind(all_mat(i).name,'_2_')
        disp = 2;
    end
    Nalpha_per_cycle = round(60/HR*1000/6);
    
    SI_tail = SI_with_motion(:,end-Nalpha_per_cycle:end,:);
    SI_min(:,:,i) = squeeze(min(SI_tail,[],2));
    SI_max(:,:,i) = squeeze(max(SI_tail,[],2));
    SI_mean(:,:,i) = squeeze(mean(SI_tail,2));
    
    HR_all(i) = HR;
    disp_all(i) = disp;
    
end
load('SI_without_motion.mat')
variation = (SI_max-SI_min)/2./SI_mean;
save('data.mat','variation','HR_all','disp_all','SI_max','SI_mean','SI_min','SI_no_motion')