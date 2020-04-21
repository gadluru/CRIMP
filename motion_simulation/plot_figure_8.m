ccc
%% parameters
ratio = 1;
flip_angle = 15; % degree
TR = 2; % [ms]
M0 = 1;
heart_rate = 120; % simulated heart rate [min^-1]
max_disp = 1; % simulated maximum myocardium displacement, normalizde by number of slice thickness
T1 = 10:10:3000; % simulated T1 range [ms]
%%
dic_real = zeros(12,100);

load('SliceFast.mat')
Slice_profile = Slice.Profile;

flip_angle = flip_angle*ratio;
flip_angle_slice = flip_angle.*Slice_profile;

Mz = ones(3800,1);
Mz_PD = ones(3800,1);
simu_T1 = 300;


Simu_T1 = simu_T1;
% assumed T1 for the entire slab
fprintf(sprintf('simulated T1 = %g ms\n', simu_T1))
fprintf(sprintf('heart rate = %g /min\n', heart_rate))
fprintf(sprintf('maximum displacement = %g slice thickness\n', max_disp))

Slab_T1 = ones(3800,1)*Simu_T1;

%% simulate Mz evolution with or without motion
[SI_no_motion,Mz_no_motion] = simu_steady_state(flip_angle_slice,Mz,M0,TR,Slab_T1);
[SI_with_motion,Mz_motion] = simu_steady_state_with_motion(flip_angle_slice,Mz,M0,TR,Slab_T1,max_disp,heart_rate);


T = 1:3000; % number of alpha pulses
total_time = T(end)*TR/1000;
Ncycle = total_time/60*heart_rate;
Motion = round(sin(T/T(end)*Ncycle*2*pi)*max_disp*100);
Motion = interp1(1:3000,Motion,1:0.5:3000);

Colors = colors;


%% Mz evolution without motion (not shown in the paper)
figure
imagesc(Mz_no_motion(801:3000,1:2000))
% axis image
colormap gray
hold on
plot([0,2000],[200,200],'Color',Colors(1,:),'LineWidth',2)
plot([0,2000],[400,400],'Color',Colors(1,:),'LineWidth',2)
plot([0,2000],[600,600],'Color',Colors(1,:),'LineWidth',2)
plot([0,2000],[800,800],'Color',Colors(1,:),'LineWidth',2)
plot([0,2000],[1000,1000],'Color',Colors(1,:),'LineWidth',2)
plot([0,2000],[1200,1200],'Color',Colors(1,:),'LineWidth',2)
plot([0,2000],[1400,1400],'Color',Colors(1,:),'LineWidth',2)
plot([0,2000],[1600,1600],'Color',Colors(1,:),'LineWidth',2)
plot([0,2000],[1800,1800],'Color',Colors(1,:),'LineWidth',2)
plot([0,2000],[2000,2000],'Color',Colors(1,:),'LineWidth',2)
yticks(300:200:1900)
yticklabels({'1','2','3','4','5','6','7','8','9'})
xlabel 'Time(s)'
ylabel 'slice number'
xticks([1,500:500:2000])
xticklabels({'0','0.5','1','1.5','2'})
set(gca,'FontSize',24)
title 'M_z evolution w/o motion'
colorbar


%% Mz evolution with motion (not shown in the paper)
figure
imagesc(Mz_motion(801:3000,1:2000))
% axis image
colormap gray
hold on
plot(Motion(1:2000)+200,'Color',Colors(1,:),'LineWidth',2)
plot(Motion(1:2000)+400,'Color',Colors(1,:),'LineWidth',2)
plot(Motion(1:2000)+600,'Color',Colors(1,:),'LineWidth',2)
plot(Motion(1:2000)+800,'Color',Colors(1,:),'LineWidth',2)
plot(Motion(1:2000)+1000,'Color',Colors(1,:),'LineWidth',2)
plot(Motion(1:2000)+1200,'Color',Colors(1,:),'LineWidth',2)
plot(Motion(1:2000)+1400,'Color',Colors(1,:),'LineWidth',2)
plot(Motion(1:2000)+1600,'Color',Colors(1,:),'LineWidth',2)
plot(Motion(1:2000)+1800,'Color',Colors(1,:),'LineWidth',2)
plot(Motion(1:2000)+2000,'Color',Colors(1,:),'LineWidth',2)
yticks(300:200:1900)
yticklabels({'1','2','3','4','5','6','7','8','9'})
xlabel 'Time(s)'
ylabel 'slice number'
xticks([1,500:500:2000])
xticklabels({'0','0.5','1','1.5','2'})
set(gca,'FontSize',24)
title 'M_z evolution w/ motion'
colorbar


%% slice profile (Supporting Information Figure S7)
figure,
plot(Slice_profile,'LineWidth',2,'Color','k')
hold on
plot([400 400],[0 1],'LineWidth',2,'Color',Colors(1,:),'LineStyle','--')
plot([600 600],[0 1],'LineWidth',2,'Color',Colors(1,:),'LineStyle','--')
plot([200 200],[0 1],'LineWidth',2,'Color',Colors(1,:),'LineStyle','--')
plot([800 800],[0 1],'LineWidth',2,'Color',Colors(1,:),'LineStyle','--')
set(gca,'FontSize',24)
xticks([])
title 'Slice Profile'
set(gcf,'Position',[1000 896 560*0.8 420*0.8])


%% signal-time curve (Figure 8a)
figure
plot(SI_no_motion(1,1:334)','LineStyle','--','LineWidth',2)
hold on
plot(SI_no_motion(2,1:334)','LineStyle','--','LineWidth',2)
plot(SI_no_motion(3,1:334)','LineStyle','--','LineWidth',2)
axis([1 334 20 55])
xticks([1 83 167 250 334])
xticklabels({'0','0.5','1','1.5','2'})
xlabel 'Time(s)'
ylabel 'Signal Intensity [a.u.]'

hold on
plot(SI_with_motion(1,1:334)','LineWidth',2,'Color',Colors(1,:))
plot(SI_with_motion(2,1:334)','LineWidth',2,'Color',Colors(2,:))
plot(SI_with_motion(3,1:334)','LineWidth',2,'Color',Colors(3,:))

patch([251,251,333,333],[23,24.5,24.5,23],0,'EdgeColor','k','FaceColor','none','LineWidth',2)
legend({'slice 1 w/o motion','slice 2 w/o motion','slice 3 w/o motion','slice 1 w/ motion','slice 2 w/ motion','slice 3 w/ motion'})
set(gca,'FontSize',24)
set(gcf,'Position',[1000 896 560*1.25 420])


%% signal-time curve zoomed (Figure 8b)
figure
plot(SI_no_motion(2,251:334)','LineStyle','--','LineWidth',2,'Color',Colors(2,:))
hold on
plot(SI_no_motion(3,251:334)','LineStyle','--','LineWidth',2,'Color',Colors(3,:))
axis([1 84 23.5 23.95])
xticks(1:83/5:84 )
xticklabels({'1.5','1.6','1.7','1.8','1.9','2'})
xlabel 'Time(s)'
ylabel 'Signal Intensity [a.u.]'

hold on
plot(SI_with_motion(2,251:334)','LineWidth',2,'Color',Colors(2,:))
plot(SI_with_motion(3,251:334)','LineWidth',2,'Color',Colors(3,:))
plot(SI_with_motion(4,251:334)','LineWidth',2,'Color',Colors(4,:))
plot(SI_with_motion(5,251:334)','LineWidth',2,'Color',Colors(5,:))
legend({'slice 2 w/o motion','slice 3-5 w/o motion','slice 2 w/ motion','slice 3 w/ motion','slice 4 w/ motion','slice 5 w/ motion'},'Position',[0.2 0.5 0.26 0.2])
set(gca,'FontSize',24)
set(gcf,'Position',[100 100 560*1.25 420])


%% Mz profile (not shown in paper)
figure
plot(Mz_no_motion(801:3000,end),'LineWidth',2)
hold on
plot(Mz_motion(801:3000,end),'LineWidth',2)

plot([200 200],[0.3 1],'LineWidth',2,'Color','k','LineStyle','--')
plot([400 400],[0.3 1],'LineWidth',2,'Color','k','LineStyle','--')
plot([600 600],[0.3 1],'LineWidth',2,'Color','k','LineStyle','--')
plot([800 800],[0.3 1],'LineWidth',2,'Color','k','LineStyle','--')
plot([1000 1000],[0.3 1],'LineWidth',2,'Color','k','LineStyle','--')
plot([1200 1200],[0.3 1],'LineWidth',2,'Color','k','LineStyle','--')
plot([1400 1400],[0.3 1],'LineWidth',2,'Color','k','LineStyle','--')
plot([1600 1600],[0.3 1],'LineWidth',2,'Color','k','LineStyle','--')
plot([1800 1800],[0.3 1],'LineWidth',2,'Color','k','LineStyle','--')
plot([2000 2000],[0.3 1],'LineWidth',2,'Color','k','LineStyle','--')
axis([1 2200 0.3 1])
ylabel 'Normalized Mz'
xticks(300:200:1900)
xticklabels({'slice 1','slice 2','slice 3','slice 4','slice 5','slice 6','slice 7','slice 8','slice 9'})
set(gca,'FontSize',12)
legend({'M_z w/o motion','M_z w/ motion'})
title('M_z profile at time=2s')
set(gcf,'Position',[1000 896 560*2 420*0.8])