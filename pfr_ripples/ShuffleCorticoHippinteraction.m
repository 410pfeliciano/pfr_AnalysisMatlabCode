Fs = 1000;
sint = 1/Fs;

%% Cortical activity
CTXLFP = ch42Tet13LFP;
ctxLFPmualfp = bandFilter (CTXLFP, 'highfreq'); % 100-450Hz Filter
ctxMUAenvp = abs(hilbert(ctxLFPmualfp));
ctxMUAenvpsm = smoothdata(ctxMUAenvp,'gaussian',20); 
%% generating random ripple indices
time1 = 1600;
time2 = 1609;
startIdx = time1/sint;
endIdx = time2/sint;
shuffleIdx = (endIdx-startIdx).*rand(10000,1) + startIdx;

%% Obtaining cortical activity
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
trials = 1:length(shuffleIdx);
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
ctxMUAmatrix = zeros(size(shuffleIdx,1),Veclenght);
for kk = 1: size(shuffleIdx,1)-1
  ctxMUAmatrix(kk,:) = ctxMUAenvpsm(shuffleIdx(kk) - eventExt : shuffleIdx(kk)+ eventExt);
end

%%

figure('DefaultAxesFontSize',10)
imagesc(tt,trials,ctxMUAmatrix)
c = colorbar;
c.Label.String = 'Ctx MUA envelope(uV)';
c.FontSize = 11;
% colorbar('southoutside')
axis xy
ylim([0 inf])
xlabel('Time relative to Shuffle Idx (sec)');
ylabel('Shuffle Event #')
caxis([20 100])
colormap hot
set(gca,'FontSize',11)
% colorbar('off')

%%
close all
% Adding the function jbfill to the path
% addpath(genpath('C:\Users\pdrfl\MATLAB Drive\pfr_AnalysisScripts\pfr_plots'));
ctxMUAavg = mean(ctxMUAmatrix,1);
ctxMUAerror = mad(ctxMUAmatrix,1);
hold on
ctx = plot(tt, ctxMUAavg, 'r', 'LineWidth',2);
plot(tt, ctxMUAavg + ctxMUAerror, '--r')
plot(tt, ctxMUAavg - ctxMUAerror, '--r')
jbfill(tt,[ctxMUAavg + ctxMUAerror],[ctxMUAavg - ctxMUAerror],... 
    'r','r',0,0.1);
line([0 0],[0,70], 'Color','k','LineStyle','--')
legend(ctx, 'Cortex MUA')
legend boxoff
xlabel('Time relative to FG peak (sec)')
ylabel('MUA envp (uV)')
axis([-inf inf 0 80])
set(gca,'FontSize',11)
hold off


%% Calculate Cortical activity before, during, and after
ctxBelowmeddian = zeros(size(shuffleIdx,1),1);
ctxCentermeddian = zeros(size(shuffleIdx,1),1);
ctxAbovmeddian = zeros(size(shuffleIdx,1),1);
for kk = 1:size(shuffleIdx,1)
    ctxBelowmeddian(kk) = median(ctxMUAmatrix(kk,1:50));
    ctxCentermeddian(kk) = median(ctxMUAmatrix(kk,75:125));
    ctxAbovmeddian(kk) = median(ctxMUAmatrix(kk,150:201));
end

figure
hold on
h= histogram(ctxBelowmeddian,'Normalization','probability',... 
    'BinWidth', 2,'BinEdges',0:100,'DisplayStyle','stairs', 'LineWidth',2);
h1= histogram(ctxCentermeddian,'Normalization','probability',... 
    'BinWidth', 2,'BinEdges',0:100,'DisplayStyle','stairs', 'LineWidth',2);
h2= histogram(ctxAbovmeddian,'Normalization','probability',... 
    'BinWidth', 2,'BinEdges',0:100,'DisplayStyle','stairs', 'LineWidth',2);
hold off
% h2.BinEdges = 0:100;
h.BinWidth  = 2;
h1.BinWidth = 2;
h2.BinWidth = 2;
% h.EdgeColor = 'k';
xlabel('Ctx MUA envp (uV)')
ylabel('Probability Density')
ylim([0 inf])
legend([h, h1, h2], 'Ctx MUA Before DG shuffleIdx',...
    'Ctx MUA During DG shuffleIdx','Ctx MUA After DG shuffleIdx')
legend('boxoff')  
set(gca,'FontSize',11)
box off

%% Classification of shuffleIdx occurence during "up-states"
earlyUp = zeros(size(shuffleIdx,1),1);
lateUp = zeros(size(shuffleIdx,1),1);
betweenUp = zeros(size(shuffleIdx,1),1);
upUp = zeros(size(shuffleIdx,1),1);

for kk = 1: size(shuffleIdx,1)
   if  ctxCentermeddian(kk) > 24 % if shuffleIdx occur during upstate
     if   ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) < 24
         betweenUp(kk) = shuffleIdx(kk);
     elseif ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) > 24
         earlyUp(kk) = shuffleIdx(kk);
     elseif ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) < 24
         lateUp(kk) = shuffleIdx(kk);
     elseif ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) > 24
         upUp(kk) = shuffleIdx(kk); 
     end
   end
end
 
earlyUp(earlyUp == 0) = [];
lateUp(lateUp == 0) = [];
betweenUp(betweenUp == 0) = [];
upUp(upUp == 0) = [];

%% Matrix calculation for shuffleIdxs occuring during Up-states
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
earlyUpMUAmat = zeros(size(earlyUp,1),Veclenght);
for kk = 1: size(earlyUp,1)
    earlyUpMUAmat(kk,:) = ctxMUAenvpsm(earlyUp(kk) - eventExt : earlyUp(kk)+ eventExt);
   
end

lateUpMUAmat= zeros(size(lateUp,1),Veclenght);
for kk = 1: size(lateUp,1)
    lateUpMUAmat(kk,:) = ctxMUAenvpsm(lateUp(kk) - eventExt : lateUp(kk)+ eventExt);
end

betweenUpMUAmat= zeros(size(betweenUp,1),Veclenght);
for kk = 1: size(betweenUp,1)
    betweenUpMUAmat(kk,:) = ctxMUAenvpsm(betweenUp(kk) - eventExt : betweenUp(kk)+ eventExt);
end

upUpMUAmat= zeros(size(upUp,1),Veclenght);
for kk = 1: size(upUp,1)
    upUpMUAmat(kk,:) = ctxMUAenvpsm(upUp(kk) - eventExt : upUp(kk)+ eventExt);
end
close all
subplot(141)
trials = 1:length(earlyUp);
imagesc(tt,trials,earlyUpMUAmat)
c = colorbar;
axis xy
ylim([0 inf])
% xlabel('Time relative to FG peak (sec)');
ylabel('shuffleIdx event #')
caxis([20 100])
colormap hot
colorbar('off')
set(gca,'FontSize',11)

subplot(143)
trials = 1:length(lateUp);
imagesc(tt,trials,lateUpMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
% xlabel('Time relative to FG peak (sec)');
% ylabel('shuffleIdx event #')
caxis([20 100])
colormap hot
colorbar('off')
set(gca,'FontSize',11)

subplot(142)
trials = 1:length(betweenUp);
imagesc(tt,trials,betweenUpMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
% xlabel('Time relative to FG peak (sec)');
% ylabel('shuffleIdx event #')
caxis([20 100])
colormap hot
colorbar('off')
set(gca,'FontSize',11)

subplot(144)
trials = 1:length(upUp);
imagesc(tt,trials,upUpMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
c.FontSize = 11;
axis xy
ylim([0 inf])
% xlabel('Time relative to Shuffle Idx');
% ylabel('shuffleIdx event #')
caxis([20 100])
colormap hot
colorbar('off')
set(gca,'FontSize',11)
%%
% lateUpMUAmatctxMUAmean = mean(lateUpMUAmat,1);
% lateUpMUAmatctxMUAerror = mad(lateUpMUAmat,1);
% close all
% hold on
% ctx = plot(tt, lateUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
% jbfill(tt,lateUpMUAmatctxMUAmean + lateUpMUAmatctxMUAerror,... 
%     lateUpMUAmatctxMUAmean - lateUpMUAmatctxMUAerror,... 
%     [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
% line([0 0],[0,100], 'Color','k','LineStyle','--')
% legend(ctx,'Ctx MUA')
% ylabel('MUA envelope (uV)')
% xlabel('Time relative to Shuffle Idx (sec)');
% legend boxoff
% set(gca,'FontSize',11)
% hold off



%
% earlyUpMUAmatctxMUAmean = mean(earlyUpMUAmat,1);
% earlyUpMUAmatctxMUAerror = mad(earlyUpMUAmat,1);
% close all
% hold on
% ctx = plot(tt, earlyUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
% jbfill(tt,earlyUpMUAmatctxMUAmean + earlyUpMUAmatctxMUAerror,... 
%     earlyUpMUAmatctxMUAmean - earlyUpMUAmatctxMUAerror,... 
%     [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
% line([0 0],[0,100], 'Color','k','LineStyle','--')
% legend(ctx,'Ctx MUA')
% ylabel('MUA envelope (uV)')
% xlabel('Time relative to Shuffle Idx (sec)');
% legend boxoff
% set(gca,'FontSize',11)
% hold off


%
% betweenUpMUAmatctxMUAmean = mean(betweenUpMUAmat,1);
% betweenUpMUAmatctxMUAerror = mad(betweenUpMUAmat,1);
% close all
% hold on
% ctx = plot(tt, betweenUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
% jbfill(tt,betweenUpMUAmatctxMUAmean + betweenUpMUAmatctxMUAerror,... 
%     betweenUpMUAmatctxMUAmean - betweenUpMUAmatctxMUAerror,... 
%     [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
% line([0 0],[0,100], 'Color','k','LineStyle','--')
% legend(ctx,'Ctx MUA')
% ylabel('MUA envelope (uV)')
% xlabel('Time relative to Shuffle Idx (sec)');
% legend boxoff
% set(gca,'FontSize',11)
% hold off


%%
% upUpMUAmatctxMUAmean = mean(upUpMUAmat,1);
% upUpMUAmatctxMUAerror = mad(upUpMUAmat,1);
% close all
% hold on
% ctx = plot(tt,upUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
% jbfill(tt,upUpMUAmatctxMUAmean +upUpMUAmatctxMUAerror,... 
%    upUpMUAmatctxMUAmean -upUpMUAmatctxMUAerror,... 
%     [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
% line([0 0],[0,100], 'Color','k','LineStyle','--')
% legend(ctx,'Ctx MUA')
% xlabel('Time relative to Shuffle Idx (sec)');
% ylabel('MUA envelope (uV)')
% legend boxoff
% set(gca,'FontSize',11)
% hold off



%% Clasification of shuffleIdx occuridng during "Down-States"
earlyDown = zeros(size(shuffleIdx,1),1);
lateDown = zeros(size(shuffleIdx,1),1);
betweenDown = zeros(size(shuffleIdx,1),1);
DownDown = zeros(size(shuffleIdx,1),1);

for kk = 1: size(shuffleIdx,1)
   if  ctxCentermeddian(kk) < 24 % if shuffleIdx occur during upstate
     if   ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) > 24
         betweenDown(kk) = shuffleIdx(kk);
     elseif ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) < 24
         earlyDown(kk) = shuffleIdx(kk);
     elseif ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) > 24
         lateDown(kk) = shuffleIdx(kk);
     elseif ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) < 24
         DownDown(kk) = shuffleIdx(kk); 
     end
   end
end
 
earlyDown(earlyDown == 0) = [];
lateDown(lateDown == 0) = [];
betweenDown(betweenDown == 0) = [];
DownDown(DownDown == 0) = [];

%% Matrix calculation for shuffleIdxs occuring during Down-states
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
earlyDownMUAmat = zeros(size(earlyDown,1),Veclenght);
for kk = 1: size(earlyDown,1)
    earlyDownMUAmat(kk,:) = ctxMUAenvpsm(earlyDown(kk) - eventExt : earlyDown(kk)+ eventExt);
end

lateDownMUAmat= zeros(size(lateDown,1),Veclenght);
for kk = 1: size(lateDown,1)
    lateDownMUAmat(kk,:) = ctxMUAenvpsm(lateDown(kk) - eventExt : lateDown(kk)+ eventExt);
end

betweenDownMUAmat= zeros(size(betweenDown,1),Veclenght);
for kk = 1: size(betweenDown,1)
    betweenDownMUAmat(kk,:) = ctxMUAenvpsm(betweenDown(kk) - eventExt : betweenDown(kk)+ eventExt);
end

DownDownMUAmat= zeros(size(DownDown,1),Veclenght);
for kk = 1: size(DownDown,1)
    DownDownMUAmat(kk,:) = ctxMUAenvpsm(DownDown(kk) - eventExt : DownDown(kk)+ eventExt);
end
% close all
% subplot(141)
% trials = 1:length(earlyDown);
% imagesc(tt,trials,earlyDownMUAmat)
% c = colorbar;
% axis xy
% ylim([0 inf])
% % xlabel('Time relative to FG peak (sec)');
% ylabel('shuffleIdx event #')
% caxis([20 100])
% colormap hot
% colorbar('off')
% set(gca,'FontSize',11)
% ylim([1 inf])
% 
% subplot(143)
% trials = 1:length(lateDown);
% imagesc(tt,trials,lateDownMUAmat)
% c = colorbar;
% c.Label.String = 'MUA spk/sec';
% axis xy
% ylim([0 inf])
% % xlabel('Time relative to FG peak (sec)');
% % ylabel('shuffleIdx event #')
% caxis([20 100])
% colormap hot
% colorbar('off')
% set(gca,'FontSize',11)
% ylim([1 inf])
% 
% subplot(142)
% trials = 1:length(betweenDown);
% imagesc(tt,trials,betweenDownMUAmat)
% c = colorbar;
% c.Label.String = 'MUA spk/sec';
% axis xy
% ylim([0 inf])
% ylim([0 inf])
% % xlabel('Time relative to FG peak (sec)');
% % ylabel('shuffleIdx event #')
% caxis([20 100])
% colormap hot
% colorbar('off')
% set(gca,'FontSize',11)
% ylim([1 inf])
% 
% subplot(144)
% trials = 1:length(DownDown);
% imagesc(tt,trials,DownDownMUAmat)
% c = colorbar;
% c.Label.String = 'MUA envelope (uV)';
% c.FontSize = 11;
% axis xy
% ylim([1 inf])
% % xlabel('Time relative to FG peak (sec)');
% % ylabel('shuffleIdx event #')
% caxis([20 100])
% ylim([1 length(trials)])
% colormap hot
% colorbar('off')
% set(gca,'FontSize',11)

%%
% lateDownMUAmatctxMUAmean = mean(lateDownMUAmat,1);
% lateDownMUAmatctxMUAerror = mad(lateDownMUAmat,1);
% close all
% hold on
% ctx = plot(tt, lateDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
% jbfill(tt,lateDownMUAmatctxMUAmean + lateDownMUAmatctxMUAerror,... 
%     lateDownMUAmatctxMUAmean - lateDownMUAmatctxMUAerror,... 
%     [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
% line([0 0],[0,100], 'Color','k','LineStyle','--')
% legend(ctx,'Ctx MUA')
% ylabel('MUA envelope (uV)')
% xlabel('Time relative to Shuffle Idx (sec)');
% legend boxoff
% set(gca,'FontSize',11)
% hold off


%
% earlyDownMUAmatctxMUAmean = mean(earlyDownMUAmat,1);
% earlyDownMUAmatctxMUAerror = mad(earlyDownMUAmat,1);
% close all
% hold on
% ctx = plot(tt, earlyDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
% jbfill(tt,earlyDownMUAmatctxMUAmean + earlyDownMUAmatctxMUAerror,... 
%     earlyDownMUAmatctxMUAmean - earlyDownMUAmatctxMUAerror,... 
%     [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
% line([0 0],[0,100], 'Color','k','LineStyle','--')
% legend(ctx,'Ctx MUA')
% ylabel('MUA envelope (uV)')
% xlabel('Time relative to Shuffle Idx (sec)');
% legend boxoff
% set(gca,'FontSize',11)
% hold off

%
% betweenDownMUAmatctxMUAmean = mean(betweenDownMUAmat,1);
% betweenDownMUAmatctxMUAerror = mad(betweenDownMUAmat,1);
% close all
% hold on
% ctx = plot(tt, betweenDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
% jbfill(tt,betweenDownMUAmatctxMUAmean + betweenDownMUAmatctxMUAerror,... 
%     betweenDownMUAmatctxMUAmean - betweenDownMUAmatctxMUAerror,... 
%     [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
% line([0 0],[0,100], 'Color','k','LineStyle','--')
% legend(ctx,'Ctx MUA')
% ylabel('MUA envelope (uV)')
% xlabel('Time relative to Shuffle Idx (sec)');
% legend boxoff
% set(gca,'FontSize',11)
% hold off


%
% DownDownMUAmatctxMUAmean = mean(DownDownMUAmat,1);
% DownDownMUAmatctxMUAerror = mad(DownDownMUAmat,1);
% close all
% hold on
% ctx = plot(tt, DownDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
% jbfill(tt,DownDownMUAmatctxMUAmean + DownDownMUAmatctxMUAerror,... 
%     DownDownMUAmatctxMUAmean - DownDownMUAmatctxMUAerror,... 
%     [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
% line([0 0],[0,100], 'Color','k','LineStyle','--')
% legend(ctx,'Ctx MUA')
% ylabel('MUA envelope (uV)')
% xlabel('Time relative to Shuffle Idx (sec)');
% legend boxoff
% set(gca,'FontSize',11)
% hold off

%%
totalRipdetected = size(lateUpMUAmat,1) + size(earlyUpMUAmat,1) ...
    + size(betweenUpMUAmat,1) + size(upUpMUAmat,1) + ...
    size(lateDownMUAmat,1) + size(earlyDownMUAmat,1) ...
    + size(betweenDownMUAmat,1) + size(DownDownMUAmat,1);
shuffleIdxUpStatePorcentage = (size(lateUpMUAmat,1) + size(earlyUpMUAmat,1) ...
    + size(betweenUpMUAmat,1) + size(upUpMUAmat,1))/totalRipdetected;
shuffleIdxDownStatePorcentage = (size(lateDownMUAmat,1) + size(earlyDownMUAmat,1) ...
    + size(betweenDownMUAmat,1) + size(DownDownMUAmat,1))/totalRipdetected;
%%
close all
X = categorical({'"Up-State"','"Down-State"'});
X = reordercats(X,{'"Up-State"','"Down-State"'});
b = bar(X,[shuffleIdxUpStatePorcentage*100, shuffleIdxDownStatePorcentage*100]);
b.FaceColor = 'flat';
b.CData(2,:) = [213/255 29/255 91/255];
b.CData(1,:) = [10/255 1/255 171/255];
ylabel('shuffleIdx Occurrence %')
ylim([0 100])
set(gca,'FontSize',11)
box off
%%
close all
subplot(121)
X = categorical({'Early Up-State','Late Up-State',...
    'Between Down-States', 'Up-States'});
X = reordercats(X,{'Early Up-State','Late Up-State',...
    'Between Down-States', 'Up-States'});
b = bar(X,[size(earlyUpMUAmat,1)/totalRipdetected*100,...
    size(lateUpMUAmat,1)/totalRipdetected*100,...
    size(betweenUpMUAmat,1)/totalRipdetected*100,...
    size(upUpMUAmat,1)/totalRipdetected*100],...
    'FaceColor', [10/255 1/255 171/255]);
b.FaceColor = 'flat';
ylabel('shuffleIdx Occurrence %')
ylim([0 70])
set(gca,'FontSize',11)
box off
subplot(122)
X = categorical({'Early Down-State','Late Down-State',...
    'Between Down-States', 'Down-States'});
X = reordercats(X,{'Early Down-State','Late Down-State',...
    'Between Down-States', 'Down-States'});
b = bar(X,[size(earlyDownMUAmat,1)/totalRipdetected*100,...
    size(lateDownMUAmat,1)/totalRipdetected*100,...
    size(betweenDownMUAmat,1)/totalRipdetected*100,...
    size(DownDownMUAmat,1)/totalRipdetected*100],...
    'FaceColor', [213/255 29/255 91/255]);
b.FaceColor = 'flat';
ylabel('shuffleIdx Occurrence %')
set(gca,'FontSize',11)
ylim([0 70])
box off











