clear; clc;
TR = 2; % [ms]

SRT_all = 0:10:200;
flip_angle = 15;

for n = 1:6
    for iSRT=1:length(SRT_all)
        
        % TR = kSpace_info.Portocol.alTR(2)/1000;
        T1 = [225,335].';
        nor = 30;
        M0 = 1;
        clear S
        S(1,1) = 0;
        S(2,1) = 0;
        SRT = SRT_all(iSRT);
        S(:,end+1) = M0 - (M0 - S(:,end)).*exp(-SRT./T1);
        clear signal
        signal(1,1) = 0;
        signal(2,1) = 0;
        for i=1:nor
            signal(:,end+1) = S(:,end).*sin(flip_angle/180*pi);
            S(:,end+1) = S(:,end).*cos(flip_angle/180*pi);
            S(:,end+1) = M0 - (M0 - S(:,end)).*exp(-(TR*n)./T1);
        end
        
        % figure,plot(S(:,3:3:end).')
        
        %Contrast(iSRT,n) = (mean(signal(1,2:end))-mean(signal(2,2:end)));
        %Signal(:,iSRT,iFA) = mean(S,2)*sin(flip_angle);
        Signal_SR(:,iSRT,n) = mean(signal,2);
    end
    
    nor = 3000;
    
    M0 = 1;
    clear S
    S(1,1) = M0;
    S(2,1) = M0;
    clear signal
    signal(1,1) = 0;
    signal(2,1) = 0;
    for i=1:nor
        signal(:,end+1) = S(:,end).*sin(flip_angle/180*pi);
        S(:,end+1) = S(:,end).*cos(flip_angle/180*pi);
        S(:,end+1) = M0 - (M0 - S(:,end)).*exp(-(TR*n)./T1);
    end
    %Contrast_SPGR(:,n) = (signal(1,end)-signal(2,end));
    Signal_SPGR(:,n) = signal(:,end);
end

color1 = [ 0    0.4470    0.7410];
color2 = [ 0.8500    0.3250    0.0980];
color3 = [0.9290    0.6940    0.1250];
color4 = [    0.4940    0.1840    0.5560];
color5 = [0.4660    0.6740    0.1880];
color6 =[0.3010    0.7450    0.9330];
color7 = [0.6350    0.0780    0.1840];

LineWidth = 3;

figure,plot(SRT_all,squeeze(Signal_SR(1,:,1:4)),'LineWidth',LineWidth)
hold on
plot(SRT_all,repmat(Signal_SPGR(1,1),[length(SRT_all),1]),'Color',color1,'LineStyle','--','LineWidth',LineWidth)
plot(SRT_all,repmat(Signal_SPGR(1,2),[length(SRT_all),1]),'Color',color2,'LineStyle','--','LineWidth',LineWidth)
plot(SRT_all,repmat(Signal_SPGR(1,3),[length(SRT_all),1]),'Color',color3,'LineStyle','--','LineWidth',LineWidth)
plot(SRT_all,repmat(Signal_SPGR(1,4),[length(SRT_all),1]),'Color',color4,'LineStyle','--','LineWidth',LineWidth)
%plot(SRT_all,repmat(Signal_SPGR(1,5),[length(SRT_all),1]),'Color',color5,'LineStyle','--','LineWidth',1.5')
%plot(SRT_all,repmat(Signal_SPGR(1,6),[length(SRT_all),1]),'Color',color6,'LineStyle','--','LineWidth',1.5')
%plot(SRT_all,repmat(Contrast_SPGR(:,5),[41,5]),'Color',color5)

legend('SR N_{SG}=1','SR N_{SG}=2','SR N_{SG}=3','SR N_{SG}=4','CRIMP N_{SG}=1','CRIMP N_{SG}=2','CRIMP N_{SG}=3','CRIMP N_{SG}=4','Location','best','NumColumns',2)
%legend('SR prepared SG=1','SR prepared SG=2','SR prepared SG=3','SR prepared SG=4','SR prepared SG=5','SR prepared SG=6','SPGR SG=1','SPGR SG=2','SPGR SG=3','SPGR SG=4','SPGR SG=5','SPGR SG=6')

grid on
set(gcf,'Position',[0,0,800,600])

xlabel 'SRT (ms)'
ylabel 'Absolute SI [a.u.]'
set(gca,'FontSize',28)
title('SI (T_1=225ms, FA=15^o)')
