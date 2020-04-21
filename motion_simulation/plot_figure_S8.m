ccc

load('data.mat')

LineWidth = 3;

[~,order] = sort(HR_all+disp_all);


for i=1:9
    d = SI_max(i,:,:) - SI_no_motion(i,:)';
    d = abs(d);
    [~,T1_max(i,:,:)] = min(d,[],1);
    
    d = SI_min(i,:,:) - SI_no_motion(i,:)';
    d = abs(d);
    [~,T1_min(i,:,:)] = min(d,[],1);
    
    d = SI_mean(i,:,:) - SI_no_motion(i,:)';
    d = abs(d);
    [~,T1_mean(i,:,:)] = min(d,[],1);
end


T1_normal = simu_T1;

Colors = colors;
N = 48; % HR=120 disp = 2
slice = 1;
%%
fig1 = figure;
T1_max_temp = interp1(10:10:3000,T1_min(slice,:,N),10:1:3000);
T1_min_temp = interp1(10:10:3000,T1_max(slice,:,N),10:1:3000);
T1_max_temp = [1:9,T1_max_temp];
T1_min_temp = [1:9,T1_min_temp];

Y = [T1_min_temp',T1_max_temp'-T1_min_temp'];
h1 = area(Y,'LineStyle','none');
h1(1).FaceAlpha = 0;
h1(2).FaceColor = (1-Colors(1,:))*0.6+Colors(1,:);

% h1(2).FaceAlpha = 1;

hold on

slice = 2;
T1_max_temp = interp1(10:10:3000,T1_min(slice,:,N),10:1:3000);
T1_min_temp = interp1(10:10:3000,T1_max(slice,:,N),10:1:3000);
T1_max_temp = [1:9,T1_max_temp];
T1_min_temp = [1:9,T1_min_temp];

Y = [T1_min_temp',T1_max_temp'-T1_min_temp'];
h2 = area(Y,'LineStyle','none');
h2(1).FaceAlpha = 0;
h2(2).FaceColor = ((1-Colors(2,:))*0.6+Colors(2,:))*0.5+((1-Colors(1,:))*0.6+Colors(1,:))*0.5;
hlegend = area(3000,1,'LineStyle','none');
hlegend.FaceColor = ((1-Colors(2,:))*0.6+Colors(2,:));
% h2(2).FaceAlpha = 1;

slice = 5;
T1_max_temp = interp1(10:10:3000,T1_min(slice,:,N),10:1:3000);
T1_min_temp = interp1(10:10:3000,T1_max(slice,:,N),10:1:3000);
T1_max_temp = [1:9,T1_max_temp];
T1_min_temp = [1:9,T1_min_temp];

Y = [T1_min_temp',T1_max_temp'-T1_min_temp'];
h3 = area(Y,'LineStyle','none');
h3(1).FaceAlpha = 0;
h3(2).FaceColor = (1-Colors(3,:))*0.6+Colors(3,:);
% h3(2).FaceAlpha = 0.4;

axis image

l0 = plot(simu_T1,simu_T1,'LineWidth',LineWidth,'Color','k');
slice = 1;
l1 = plot(T1_normal,T1_mean(slice,:,N),'Color',Colors(1,:),'LineWidth',LineWidth);
slice = 2;
l2 = plot(T1_normal,T1_mean(slice,:,N),'Color',Colors(2,:),'LineWidth',LineWidth);
slice = 5;
l3 = plot(T1_normal,T1_mean(slice,:,N),'Color',Colors(3,:),'LineWidth',LineWidth);


axis image
axis([1 1500 1 1800])

xticks([1,300,600,900,1200,1500])
yticks([1,300,600,900,1200,1500,1800])

legend([l0,l1,h1(2),l2,hlegend,l3,h3(2)],'ideal condition','slice 1','slice 1 error','slice 2','slice 2 error','slice 5','slice 5 error','Location','northwest')
ylabel 'Fitted T_1 (ms)'
xlabel 'Simulated T_1 (ms)'
title 'HR=120/min, Max. Disp.=2\times[ST]'
set(gcf,'Position',[400,600,560,560])
set(gca,'FontSize',20)
grid on
set(gca,'Layer','top')
hgexport(fig1,'FIt_error_HR120_Disp2.eps');
%%


Colors = colors;
N = 95;% HR=80 disp = 1
slice = 1;

fig2 = figure;
T1_max_temp = interp1(10:10:3000,T1_min(slice,:,N),10:1:3000);
T1_min_temp = interp1(10:10:3000,T1_max(slice,:,N),10:1:3000);
T1_max_temp = [1:9,T1_max_temp];
T1_min_temp = [1:9,T1_min_temp];

Y = [T1_min_temp',T1_max_temp'-T1_min_temp'];
h1 = area(Y,'LineStyle','none');
h1(1).FaceAlpha = 0;
h1(2).FaceColor = (1-Colors(1,:))*0.6+Colors(1,:);

% h1(2).FaceAlpha = 1;

hold on

slice = 2;
T1_max_temp = interp1(10:10:3000,T1_min(slice,:,N),10:1:3000);
T1_min_temp = interp1(10:10:3000,T1_max(slice,:,N),10:1:3000);
T1_max_temp = [1:9,T1_max_temp];
T1_min_temp = [1:9,T1_min_temp];

Y = [T1_min_temp',T1_max_temp'-T1_min_temp'];
h2 = area(Y,'LineStyle','none');
h2(1).FaceAlpha = 0;
h2(2).FaceColor = ((1-Colors(2,:))*0.6+Colors(2,:))*0.5+((1-Colors(1,:))*0.6+Colors(1,:))*0.5;
hlegend = area(3000,1,'LineStyle','none');
hlegend.FaceColor = ((1-Colors(2,:))*0.6+Colors(2,:));
% h2(2).FaceAlpha = 1;


slice = 5;
T1_max_temp = interp1(10:10:3000,T1_min(slice,:,N),10:1:3000);
T1_min_temp = interp1(10:10:3000,T1_max(slice,:,N),10:1:3000);
T1_max_temp = [1:9,T1_max_temp];
T1_min_temp = [1:9,T1_min_temp];

Y = [T1_min_temp',T1_max_temp'-T1_min_temp'];
h3 = area(Y,'LineStyle','none');
h3(1).FaceAlpha = 0;
h3(2).FaceColor = (1-Colors(3,:))*0.6+Colors(3,:);
% h3(2).FaceAlpha = 0.4;

axis image

l0 = plot(simu_T1,simu_T1,'LineWidth',LineWidth,'Color','k');
slice = 1;
l1 = plot(T1_normal,T1_mean(slice,:,N),'Color',Colors(1,:),'LineWidth',LineWidth);
slice = 2;
l2 = plot(T1_normal,T1_mean(slice,:,N),'Color',Colors(2,:),'LineWidth',LineWidth);
slice = 5;
l3 = plot(T1_normal,T1_mean(slice,:,N),'Color',Colors(3,:),'LineWidth',LineWidth);


axis image
axis([1 1500 1 1800])

xticks([1,300,600,900,1200,1500])
yticks([1,300,600,900,1200,1500,1800])

legend([l0,l1,h1(2),l2,hlegend,l3,h3(2)],'ideal condition','slice 1','slice 1 error','slice 2','slice 2 error','slice 5','slice 5 error','Location','northwest')
ylabel 'Fitted T_1 (ms)'
xlabel 'Simulated T_1 (ms)'
title 'HR=80/min, Max. Disp.=1\times[ST]'
set(gcf,'Position',[400,600,560,560])
set(gca,'FontSize',20)
grid on
set(gca,'Layer','top')
hgexport(fig2,'FIt_error_HR80_Disp1.eps');

%%

Colors = colors;
N = 49; % HR=60 disp = 0.5
slice = 1;

fig3 = figure;
T1_max_temp = interp1(10:10:3000,T1_min(slice,:,N),10:1:3000);
T1_min_temp = interp1(10:10:3000,T1_max(slice,:,N),10:1:3000);
T1_max_temp = [1:9,T1_max_temp];
T1_min_temp = [1:9,T1_min_temp];

Y = [T1_min_temp',T1_max_temp'-T1_min_temp'];
h1 = area(Y,'LineStyle','none');
h1(1).FaceAlpha = 0;
h1(2).FaceColor = (1-Colors(1,:))*0.6+Colors(1,:);

% h1(2).FaceAlpha = 1;

hold on

slice = 2;
T1_max_temp = interp1(10:10:3000,T1_min(slice,:,N),10:1:3000);
T1_min_temp = interp1(10:10:3000,T1_max(slice,:,N),10:1:3000);
T1_max_temp = [1:9,T1_max_temp];
T1_min_temp = [1:9,T1_min_temp];

Y = [T1_min_temp',T1_max_temp'-T1_min_temp'];
h2 = area(Y,'LineStyle','none');
h2(1).FaceAlpha = 0;
h2(2).FaceColor = ((1-Colors(2,:))*0.6+Colors(2,:))*0.5+((1-Colors(1,:))*0.6+Colors(1,:))*0.5;
hlegend = area(3000,1,'LineStyle','none');
hlegend.FaceColor = ((1-Colors(2,:))*0.6+Colors(2,:));
% h2(2).FaceAlpha = 1;


slice = 5;
T1_max_temp = interp1(10:10:3000,T1_min(slice,:,N),10:1:3000);
T1_min_temp = interp1(10:10:3000,T1_max(slice,:,N),10:1:3000);
T1_max_temp = [1:9,T1_max_temp];
T1_min_temp = [1:9,T1_min_temp];

Y = [T1_min_temp',T1_max_temp'-T1_min_temp'];
h3 = area(Y,'LineStyle','none');
h3(1).FaceAlpha = 0;
h3(2).FaceColor = (1-Colors(3,:))*0.6+Colors(3,:);
% h3(2).FaceAlpha = 0.4;

axis image

l0 = plot(simu_T1,simu_T1,'LineWidth',LineWidth,'Color','k');
slice = 1;
l1 = plot(T1_normal,T1_mean(slice,:,N),'Color',Colors(1,:),'LineWidth',LineWidth);
slice = 2;
l2 = plot(T1_normal,T1_mean(slice,:,N),'Color',Colors(2,:),'LineWidth',LineWidth);
slice = 5;
l3 = plot(T1_normal,T1_mean(slice,:,N),'Color',Colors(3,:),'LineWidth',LineWidth);


axis image
axis([1 1500 1 1800])

xticks([1,300,600,900,1200,1500])
yticks([1,300,600,900,1200,1500,1800])

legend([l0,l1,h1(2),l2,hlegend,l3,h3(2)],'ideal condition','slice 1','slice 1 error','slice 2','slice 2 error','slice 5','slice 5 error','Location','northwest')
ylabel 'Fitted T_1 (ms)'
xlabel 'Simulated T_1 (ms)'
title 'HR=60/min, Max. Disp.=0.5\times[ST]'
set(gcf,'Position',[400,600,560,560])
set(gca,'FontSize',20)
grid on
set(gca,'Layer','top')
hgexport(fig3,'FIt_error_HR60_Disp0.5.eps');


%%
fig4 = figure;
alphas = 0.1:0.9/6:1;
for i=1:7
    lines = ((i-1)*16+1):i*16;

    hold on
    l1 = plot([0 1],[i,i],'Color',Colors(1,:),'LineWidth',4);
    l1.Color(4)=alphas(i);
    
    l2 = plot([0 1],[i,i]-10,'Color',Colors(2,:),'LineWidth',4);
    l2.Color(4)=alphas(i);
    
    l3 = plot([0 1],[i,i]-20,'Color',Colors(3,:),'LineWidth',4);
    l3.Color(4)=alphas(i);
    
    l4 = plot([0 1],[i,i]-30,'Color',Colors(4,:),'LineWidth',4);
    l4.Color(4)=alphas(i);
    
    l5 = plot([0 1],[i,i]-40,'Color',Colors(5,:),'LineWidth',4);
    l5.Color(4)=alphas(i);
    
    l6 = plot([0 1],[i,i]-50,'Color',Colors(6,:),'LineWidth',4);
    l6.Color(4)=alphas(i);
end
axis([-0.2 3.4 -53 10])
box on
set(gca,'XTick',[],'YTick',[])

text(1.2,4,'T1=100ms, HR=60-120/min','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',20)
text(1.2,4-10,'T1=200ms, HR=60-120/min','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',20)
text(1.2,4-20,'T1=300ms, HR=60-120/min','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',20)
text(1.2,4-30,'T1=500ms, HR=60-120/min','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',20)
text(1.2,4-40,'T1=1000ms, HR=60-120/min','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',20)
text(1.2,4-50,'T1=2000ms, HR=60-120/min','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',20)
set(gca,'FontSize',20)
hgexport(fig4,'e_Legend.eps')
%%
alphas = 0.1:0.9/6:1;
fig5 = figure;
for i=1:7
    lines = ((i-1)*16+1):i*16;

    hold on
    l1 = plot(squeeze(variation(1,10,order(lines))),'Color',Colors(1,:),'LineWidth',LineWidth);
    l1.Color(4)=alphas(i);
    
    l2 = plot(squeeze(variation(1,20,order(lines))),'Color',Colors(2,:),'LineWidth',LineWidth);
    l2.Color(4)=alphas(i);
    
    l3 = plot(squeeze(variation(1,30,order(lines))),'Color',Colors(3,:),'LineWidth',LineWidth);
    l3.Color(4)=alphas(i);
    
    l4 = plot(squeeze(variation(1,50,order(lines))),'Color',Colors(4,:),'LineWidth',LineWidth);
    l4.Color(4)=alphas(i);
    
    l5 = plot(squeeze(variation(1,100,order(lines))),'Color',Colors(5,:),'LineWidth',LineWidth);
    l5.Color(4)=alphas(i);
    
    l6 = plot(squeeze(variation(1,200,order(lines))),'Color',Colors(6,:),'LineWidth',LineWidth);
    l6.Color(4)=alphas(i);
end
% legend([l1,l2,l3,l4,l5,l6],{'T1=100ms','T1=200ms','T1=300ms','T1=500ms','T1=1000ms','T1=2000ms'},'Location','northwest');
axis([1 16 0 0.7])

xticks(2:2:16)
xticklabels(0.6:0.2:2)
yticks(0:0.1:0.7)
yticklabels({'0','10%','20%','30%','40%','50%','60%','70%'})

xlabel 'Max. Disp. [ST]'
ylabel 'Signal Variation'
title 'Slice 1'
set(gca,'FontSize',20)
box on
grid on
hgexport(fig5,'e_sl1.eps')


%%
alphas = 0.1:0.9/6:1;
fig6 = figure;
for i=1:7
    lines = ((i-1)*16+1):i*16;

    hold on
    l1 = plot(squeeze(variation(2,10,order(lines))),'Color',Colors(1,:),'LineWidth',LineWidth);
    l1.Color(4)=alphas(i);
    
    l2 = plot(squeeze(variation(2,20,order(lines))),'Color',Colors(2,:),'LineWidth',LineWidth);
    l2.Color(4)=alphas(i);
    
    l3 = plot(squeeze(variation(2,30,order(lines))),'Color',Colors(3,:),'LineWidth',LineWidth);
    l3.Color(4)=alphas(i);
    
    l4 = plot(squeeze(variation(2,50,order(lines))),'Color',Colors(4,:),'LineWidth',LineWidth);
    l4.Color(4)=alphas(i);
    
    l5 = plot(squeeze(variation(2,100,order(lines))),'Color',Colors(5,:),'LineWidth',LineWidth);
    l5.Color(4)=alphas(i);
    
    l6 = plot(squeeze(variation(2,200,order(lines))),'Color',Colors(6,:),'LineWidth',LineWidth);
    l6.Color(4)=alphas(i);
end
% lgd = legend([l1,l2,l3,l4,l5,l6],{'T1=100ms','T1=200ms','T1=300ms','T1=500ms','T1=1000ms','T1=2000ms'},'Location','northwest');
axis([1 16 0 0.2])

xticks(2:2:16)
xticklabels(0.6:0.2:2)
yticks(0:0.05:0.2)
yticklabels({'0','5%','10%','15%','20%'})

xlabel 'Max. Disp. [ST]'
ylabel 'Signal Variation'
title 'Slice 2'
set(gca,'FontSize',20)
box on
grid on
hgexport(fig6,'e_sl2.eps')

%%

alphas = 0.1:0.9/6:1;
fig7 = figure;
for i=1:7
    lines = ((i-1)*16+1):i*16;

    hold on
    l1 = plot(squeeze(variation(3,10,order(lines))),'Color',Colors(1,:),'LineWidth',LineWidth);
    l1.Color(4)=alphas(i);
    
    l2 = plot(squeeze(variation(3,20,order(lines))),'Color',Colors(2,:),'LineWidth',LineWidth);
    l2.Color(4)=alphas(i);
    
    l3 = plot(squeeze(variation(3,30,order(lines))),'Color',Colors(3,:),'LineWidth',LineWidth);
    l3.Color(4)=alphas(i);
    
    l4 = plot(squeeze(variation(3,50,order(lines))),'Color',Colors(4,:),'LineWidth',LineWidth);
    l4.Color(4)=alphas(i);
    
    l5 = plot(squeeze(variation(3,100,order(lines))),'Color',Colors(5,:),'LineWidth',LineWidth);
    l5.Color(4)=alphas(i);
    
    l6 = plot(squeeze(variation(3,200,order(lines))),'Color',Colors(6,:),'LineWidth',LineWidth);
    l6.Color(4)=alphas(i);
end
% lgd = legend([l1,l2,l3,l4,l5,l6],{'T1=100ms','T1=200ms','T1=300ms','T1=500ms','T1=1000ms','T1=2000ms'},'Location','northwest');
axis([1 16 0 0.01])
set(gca,'FontSize',20)
xticks(2:2:16)
xticklabels(0.6:0.2:2)
yticks(0:0.002:0.01)
yticklabels({'0','0.02%','0.04%','0.06%','0.08%','0.1%'})
xlabel 'Max. Disp. [ST]'
ylabel 'Signal Variation'
title 'Slice 3'
box on
grid on
hgexport(fig7,'e_sl3.eps')


%%
alphas = 0.1:0.9/6:1;
fig8 = figure;
for i=1:7
    lines = ((i-1)*16+1):i*16;

    hold on
    l1 = plot(squeeze(variation(4,10,order(lines))),'Color',Colors(1,:),'LineWidth',LineWidth);
    l1.Color(4)=alphas(i);
    
    l2 = plot(squeeze(variation(4,20,order(lines))),'Color',Colors(2,:),'LineWidth',LineWidth);
    l2.Color(4)=alphas(i);
    
    l3 = plot(squeeze(variation(4,30,order(lines))),'Color',Colors(3,:),'LineWidth',LineWidth);
    l3.Color(4)=alphas(i);
    
    l4 = plot(squeeze(variation(4,50,order(lines))),'Color',Colors(4,:),'LineWidth',LineWidth);
    l4.Color(4)=alphas(i);
    
    l5 = plot(squeeze(variation(4,100,order(lines))),'Color',Colors(5,:),'LineWidth',LineWidth);
    l5.Color(4)=alphas(i);
    
    l6 = plot(squeeze(variation(4,200,order(lines))),'Color',Colors(6,:),'LineWidth',LineWidth);
    l6.Color(4)=alphas(i);
end
% lgd = legend([l1,l2,l3,l4,l5,l6],{'T1=100ms','T1=200ms','T1=300ms','T1=500ms','T1=1000ms','T1=2000ms'},'Location','northeast');
axis([1 16 0 0.007])
set(gca,'FontSize',20)
xticks(2:2:16)
xticklabels(0.6:0.2:2)
yticks(0:0.001:0.007)
yticklabels({'0','0.01%','0.02%','0.03%','0.04%','0.05%','0.06%','0.07%'})
xlabel 'Max. Disp. [ST]'
ylabel 'Signal Variation'
title 'Slice 4'
box on
grid on
hgexport(fig8,'e_sl4.eps')
%%

alphas = 0.1:0.9/6:1;
fig9 = figure;
for i=1:7
    lines = ((i-1)*16+1):i*16;

    hold on
    l1 = plot(squeeze(variation(5,10,order(lines))),'Color',Colors(1,:),'LineWidth',LineWidth);
    l1.Color(4)=alphas(i);
    
    l2 = plot(squeeze(variation(5,20,order(lines))),'Color',Colors(2,:),'LineWidth',LineWidth);
    l2.Color(4)=alphas(i);
    
    l3 = plot(squeeze(variation(5,30,order(lines))),'Color',Colors(3,:),'LineWidth',LineWidth);
    l3.Color(4)=alphas(i);
    
    l4 = plot(squeeze(variation(5,50,order(lines))),'Color',Colors(4,:),'LineWidth',LineWidth);
    l4.Color(4)=alphas(i);
    
    l5 = plot(squeeze(variation(5,100,order(lines))),'Color',Colors(5,:),'LineWidth',LineWidth);
    l5.Color(4)=alphas(i);
    
    l6 = plot(squeeze(variation(5,200,order(lines))),'Color',Colors(6,:),'LineWidth',LineWidth);
    l6.Color(4)=alphas(i);
end
% lgd = legend([l1,l2,l3,l4,l5,l6],{'T1=100ms','T1=200ms','T1=300ms','T1=500ms','T1=1000ms','T1=2000ms'},'Location','northeast');
axis([1 16 0 0.008])
set(gca,'FontSize',20)
xticks(2:2:16)
xticklabels(0.6:0.2:2)
yticks(0:0.001:0.008)
yticklabels({'0','0.01%','0.02%','0.03%','0.04%','0.05%','0.06%','0.07%','0.08%'})
xlabel 'Max. Disp. [ST]'
ylabel 'Signal Variation'
title 'Slice 5'
box on
grid on
hgexport(fig9,'e_sl5.eps')

