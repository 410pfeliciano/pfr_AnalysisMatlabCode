% Ripple visualization
Fs = 1000;
sint = 1/Fs;
%% Visual Inspection of LFP and Ripple Candidates
% close all
time1= 4808;
time2= 5000;
xt1= time1/sint;
xt2= time2/sint;
interBelowIdxctx = find(belowDS.*sint >= time1-5 & belowDS.*sint <= time2+5);
% figure
subplot(811) % Cortical LFP with Down state 
plot(time(xt1:xt2), CTXLFP(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdxctx)
    plot(time(belowDS(interBelowIdxctx(kk)):abovDS(interBelowIdxctx(kk))),...
        CTXLFP(belowDS(interBelowIdxctx(kk)):abovDS(interBelowIdxctx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -inf inf])

subplot(812) % Cortical MUA
hold on
plot(spkClust(CtxTetIdx(1)).spkTime, ... 
    ones(length(spkClust(CtxTetIdx(1)).spkTime),1),'r.','MarkerSize',5)
for kk = 2: size(CtxTetIdx,1)
plot(spkClust(CtxTetIdx(kk)).spkTime, ...
    ones(length(spkClust(CtxTetIdx(kk)).spkTime),1) * kk,... 
    'r.','MarkerSize',5)
end
hold off
xlim([time1 time2])
ylim([0 inf])
% set(gca,'xtickLabel',[])
xlabel('Time[sec]')
ylabel('MUA #')

subplot(813) % Cortical MUA spk density
plot(time(xt1:xt2),CtxspkDensity(xt1:xt2),'r','LineWidth',2) 
xlabel('Time[sec]')
ylabel('spks/s')
axis tight
%
% subplot(814) % Cortical spectrogram
% interval = Fs;
% overlap = round(Fs*0.95);
% nfft = length(CTXLFP(xt1:xt2));
% LFP = CTXLFP(xt1:xt2);
% [S, F, T,P] = spectrogram(LFP,interval, overlap, nfft, Fs);
% imagesc(T,F,(P))
% colorbar
% % colorbar('southoutside')
% axis xy
% ylim([0 250])
% xlabel('Time [s]');
% ylabel('Frequency [Hz]')
% caxis([0 50])
% colormap hot
% colorbar('off')

% Hippocampus plotting
CA3LFP =hippLFP;
subplot(815)
interBelowIdx = find(belowIdx.*sint >= time1-5 & belowIdx.*sint <= time2+5);
plot(time(xt1:xt2), CA3LFP(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(belowIdx(interBelowIdx(kk)):abovIdx(interBelowIdx(kk))),...
        CA3LFP(belowIdx(interBelowIdx(kk)):abovIdx(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
axis([time1 time2 -inf inf])
% line([time1 time2],[th th],'Color','b','LineStyle','--','LineWidth',2)
hold off

subplot(816)
hold on
plot(spkClust(CA3TetIdx(1)).spkTime, ... 
    ones(length(spkClust(CA3TetIdx(1)).spkTime),1),'r.','MarkerSize',5)
for kk = 2: size(CA3TetIdx,1)
plot(spkClust(CA3TetIdx(kk)).spkTime, ...
    ones(length(spkClust(CA3TetIdx(kk)).spkTime),1) * kk,... 
    'r.','MarkerSize',5)
end
hold off
xlim([time1 time2])
ylim([0 inf])
% set(gca,'xtickLabel',[])
xlabel('Time[sec]')
ylabel('MUA #')

subplot(817)
plot(time(xt1:xt2),CA3spkDensity(xt1:xt2),'r','LineWidth',2) 
xlabel('Time[sec]')
ylabel('spks/s')
axis tight

% subplot(818) % Cortical spectrogram
% % interval = Fs;
% % overlap = round(Fs*0.95);
% % nfft = length(CA3LFP(xt1:xt2));
% % [S, F, T,P] = spectrogram(CA3LFP(xt1:xt2),interval, overlap, nfft, Fs);
% % imagesc(T,F,10*log10(P))
% % colorbar
% % colorbar('southoutside')
% axis xy
% ylim([60 250])
% xlabel('Time [s]');
% ylabel('Frequency [Hz]')
% caxis([0 20])
% colormap hot
% colorbar('off')

%%
% Cortical activity
smWin = 0.05; % smoothding window in seconds
ctxMUAenvp = abs(hilbert(ch42Tet13LFPhighfreq));
ctxMUAenvpsm = smoothdata(ctxMUAenvp,'gaussian',20); 
% CA3 activity
hippRippleenv = smoothdata(abs(hilbert(ch9Tet7LFPripple)),'gaussian',smWin*Fs);
ZhippRippleenv = zscore(hippRippleenv);
hippRipple = ch9Tet7LFPripple;

%% Visual Inspection of Cortical activity based on ripple detection
close all
eventExt = 0.1 / sint;
for kk = 100:120 %: size(,1)
   figure
   subplot(611)
   plot(time(belowIdx(kk) - eventExt : abovIdx(kk)+ eventExt), ...
       CTXLFP(belowIdx(kk) - eventExt : abovIdx(kk)+ eventExt),'k')
   axis([-inf inf, -inf inf])
   
   subplot(612)
   plot(time(belowIdx(kk) - eventExt : abovIdx(kk)+ eventExt),...
       ctxMUAenvpsm(belowIdx(kk) - eventExt : abovIdx(kk)+ eventExt),'k')
   axis([-inf inf, 0 inf])
   
   subplot(613)
   plot(time(belowIdx(kk) - eventExt : abovIdx(kk)+ eventExt), ...
       CA3LFP(belowIdx(kk) - eventExt : abovIdx(kk)+ eventExt),'r')
   axis([-inf inf, -inf inf])
   
   subplot(614)
   plot(time(belowIdx(kk) - eventExt : abovIdx(kk)+ eventExt), ...
       hippRipple(belowIdx(kk) - eventExt : abovIdx(kk)+ eventExt),'r')
   axis([-inf inf, -inf inf])
   
   subplot(615)
   plot(time(belowIdx(kk) - eventExt : abovIdx(kk)+ eventExt),...
       ZhippRippleenv(belowIdx(kk) - eventExt : abovIdx(kk)+ eventExt),'r')
   axis([-inf inf, -inf inf])
   
   subplot(616)
    plot(time(belowIdx(kk) - eventExt : abovIdx(kk)+ eventExt),... 
        CA3spkDensity(belowIdx(kk) - eventExt : abovIdx(kk)+ eventExt),... 
        'r','LineWidth',2) 
    xlabel('Time[sec]')
    ylabel('spks/s')
    axis tight
   axis([-inf inf, 0 inf])
end


%% Obtaining the Ripple peak
peakRipple = zeros(size(belowIdx,1),1);
for kk = 1: size(belowIdx,1)-1
   [R, peakRipple(kk)] = max(ZhippRippleenv(belowIdx(kk): abovIdx(kk)));
end
% this is to obtain the center indice based on the max ripple power peak
peakRippleIdx = zeros(size(peakRipple,1),1);
for kk = 1: size(belowIdx,1)-1
   v = belowIdx(kk): abovIdx(kk);
   peakRippleIdx(kk) = v(peakRipple(kk));
end
%% Before and after distributions
maxRipplepks = zeros(size(belowIdx,1),1);
maxMUApks = zeros(size(belowIdx,1),1);
maxMUActx = zeros(size(belowIdx,1),1);
for kk = 1: size(belowIdx,1)-1
   maxRipplepks(kk) = ZhippRippleenv(peakRippleIdx(kk));
   maxMUApks(kk) = CA3spkDensity(peakRippleIdx(kk));
   maxMUActx(kk) = ctxMUAenvpsm(peakRippleIdx(kk));
end
%% Distribution of MUA 
close all
figure
h= histogram(maxRipplepks,'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',1);
h.BinEdges = 0:15;
h.BinWidth = 1;
h.EdgeColor = 'k';
xlabel('Z-Score Ripple Power')
ylabel('Probability Density')
legend('Ripple Power')
legend('boxoff')  
box off
hold on
h= histogram(maxMUApks,'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',1);
h.BinEdges = 0:100;
h.BinWidth = 1;
h.EdgeColor = 'r';
xlabel('Max MUA ripple activation (Spk/sec)')
ylabel('Probability Density')
legend('MUA ')
legend('boxoff')  
box off
hold off
%% Eliminating ripple candidates with a standar deaviation < 4

pksRippleIdx = peakRippleIdx;
pksRippleIdx(maxRipplepks < 3) = []; % | maxMUActx > 100
ctxMUA = maxMUActx;
ctxMUA(maxRipplepks < 3 ) = []; % | maxMUActx > 100

%%
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
trials = 1:length(ctxMUA);
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
ctxMUAmatrix = zeros(size(ctxMUA,1),Veclenght);
for kk = 1: size(ctxMUA,1)
  ctxMUAmatrix(kk,:) = ctxMUAenvpsm(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt);
end
% Cortex spiking activity
ctxSPKmatrix = zeros(size(ctxMUA,1),Veclenght);
for kk = 1: size(ctxMUA,1)
  ctxSPKmatrix(kk,:) = CtxspkDensity(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt);
end

% Hippocampal MUA activity
hippMUAmatrix = zeros(size(pksRippleIdx,1),Veclenght);
for kk = 1: size(pksRippleIdx,1)
  hippMUAmatrix(kk,:) = CA3spkDensity(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt);
end
%%
% B = sort(hippMUAmatrix,1);
subplot(121)
imagesc(tt,trials,ctxMUAmatrix)
c = colorbar;
% c.Label.String = 'MUA envp(uV)';
% colorbar('southoutside')
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
ylabel('Ripple event #')
caxis([20 100])
colormap hot
% colorbar('off')
subplot(122)
imagesc(tt,trials,hippMUAmatrix)
c = colorbar;
c.Label.String = 'MUA spk/sec';
% colorbar('southoutside')
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 70])
colormap hot
% colorbar('off')


%%
close all
% Adding the function jbfill to the path
% addpath(genpath('C:\Users\pdrfl\MATLAB Drive\pfr_AnalysisScripts\pfr_plots'));
ctxMUAavg = mean(ctxMUAmatrix,1);
ctxMUAerror = mad(ctxMUAmatrix,1);
hippMUAavg = mean(hippMUAmatrix,1);
hippMUAerror = mad(hippMUAmatrix,1);
subplot(211)
hold on
ctx = plot(tt, ctxMUAavg, 'r', 'LineWidth',2);
plot(tt, ctxMUAavg + ctxMUAerror, '--r')
plot(tt, ctxMUAavg - ctxMUAerror, '--r')
jbfill(tt,[ctxMUAavg + ctxMUAerror],[ctxMUAavg - ctxMUAerror],... 
    'r','r',0,0.1);
line([0 0],[0,70], 'Color','k','LineStyle','--')
legend(ctx, 'Cortex MUA')
legend boxoff
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA envp (uV)')
axis([-inf inf 0 70])
hold off
subplot(212)
hold on
hipp = plot(tt, hippMUAavg, 'b', 'LineWidth',2);
plot(tt, hippMUAavg + hippMUAerror, '--b')
plot(tt, hippMUAavg - hippMUAerror, '--b')
jbfill(tt,hippMUAavg + hippMUAerror,hippMUAavg - hippMUAerror,... 
    'b','b',0,0.1);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA (spk/sec)')
legend(hipp, 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])

figure
hold on
ctx = plot(tt, ctxMUAavg, 'r', 'LineWidth',2);
plot(tt, ctxMUAavg + ctxMUAerror, '--r')
plot(tt, ctxMUAavg - ctxMUAerror, '--r')
jbfill(tt,[ctxMUAavg + ctxMUAerror],[ctxMUAavg - ctxMUAerror],... 
    'r','r',0,0.1);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hipp = plot(tt, hippMUAavg, 'b', 'LineWidth',2);
plot(tt, hippMUAavg + hippMUAerror, '--b')
plot(tt, hippMUAavg - hippMUAerror, '--b')
jbfill(tt,hippMUAavg + hippMUAerror,hippMUAavg - hippMUAerror,... 
    'b','b',0,0.1);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])

%% Calculate Cortical activity before, during, and after
ctxBelowmeddian = zeros(size(pksRippleIdx,1),1);
ctxCentermeddian = zeros(size(pksRippleIdx,1),1);
ctxAbovmeddian = zeros(size(pksRippleIdx,1),1);
for kk = 1:size(pksRippleIdx,1)
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
ylim([0 0.08])
xlabel('Ctx MUA envp (uV)')
ylabel('Probability Density')
legend([h, h1, h2], 'Before Ripple Ctx MUA','During Ripple Ctx MUA','After Ripple Ctx MUA')
legend('boxoff')  
box off

%% Classification of ripple occurence during "up-states"
earlyUp = zeros(size(pksRippleIdx,1),1);
lateUp = zeros(size(pksRippleIdx,1),1);
betweenUp = zeros(size(pksRippleIdx,1),1);
upUp = zeros(size(pksRippleIdx,1),1);

for kk = 1: size(pksRippleIdx,1)
   if  ctxCentermeddian(kk) > 24 % if ripple occur during upstate
     if   ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) < 24
         betweenUp(kk) = pksRippleIdx(kk);
     elseif ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) > 24
         earlyUp(kk) = pksRippleIdx(kk);
     elseif ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) < 24
         lateUp(kk) = pksRippleIdx(kk);
     elseif ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) > 24
         upUp(kk) = pksRippleIdx(kk); 
     end
   end
end
 
earlyUp(earlyUp == 0) = [];
lateUp(lateUp == 0) = [];
betweenUp(betweenUp == 0) = [];
upUp(upUp == 0) = [];

%% Matrix calculation for Ripples occuring during Up-states
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
earlyUpMUAmat = zeros(size(earlyUp,1),Veclenght);
hippearlyUpMUAmat = zeros(size(earlyUp,1),Veclenght);
hippearlyUpRipplemat = zeros(size(earlyUp,1),Veclenght);
for kk = 1: size(earlyUp,1)
    earlyUpMUAmat(kk,:) = ctxMUAenvpsm(earlyUp(kk) - eventExt : earlyUp(kk)+ eventExt);
    hippearlyUpMUAmat(kk,:) = CA3spkDensity(earlyUp(kk) - eventExt : earlyUp(kk)+ eventExt);
    hippearlyUpRipplemat(kk,:) = ZhippRippleenv(earlyUp(kk) - eventExt : earlyUp(kk)+ eventExt);
end

lateUpMUAmat= zeros(size(lateUp,1),Veclenght);
hipplateUpMUAmat = zeros(size(lateUp,1),Veclenght);
hipplateUpRipplemat = zeros(size(lateUp,1),Veclenght);
for kk = 1: size(lateUp,1)
    lateUpMUAmat(kk,:) = ctxMUAenvpsm(lateUp(kk) - eventExt : lateUp(kk)+ eventExt);
    hipplateUpMUAmat(kk,:) = CA3spkDensity(lateUp(kk) - eventExt : lateUp(kk)+ eventExt);
    hipplateUpRipplemat(kk,:) = ZhippRippleenv(lateUp(kk) - eventExt : lateUp(kk)+ eventExt);
end

betweenUpMUAmat= zeros(size(betweenUp,1),Veclenght);
hippbetweenUpMUAmat = zeros(size(betweenUp,1),Veclenght);
hippbetweenUpRipplemat = zeros(size(betweenUp,1),Veclenght);
for kk = 1: size(betweenUp,1)
    betweenUpMUAmat(kk,:) = ctxMUAenvpsm(betweenUp(kk) - eventExt : betweenUp(kk)+ eventExt);
    hippbetweenUpMUAmat(kk,:) = CA3spkDensity(betweenUp(kk) - eventExt : betweenUp(kk)+ eventExt);
    hippbetweenUpRipplemat(kk,:) = ZhippRippleenv(betweenUp(kk) - eventExt : betweenUp(kk)+ eventExt);
end

upUpMUAmat= zeros(size(upUp,1),Veclenght);
hippupUpMUAmat = zeros(size(upUp,1),Veclenght);
hippupUpRipplemat = zeros(size(upUp,1),Veclenght);
for kk = 1: size(upUp,1)
    upUpMUAmat(kk,:) = ctxMUAenvpsm(upUp(kk) - eventExt : upUp(kk)+ eventExt);
    hippupUpMUAmat(kk,:) = CA3spkDensity(upUp(kk) - eventExt : upUp(kk)+ eventExt);
    hippupUpRipplemat(kk,:) = ZhippRippleenv(upUp(kk) - eventExt : upUp(kk)+ eventExt);
end
subplot(141)
trials = 1:length(earlyUp);
imagesc(tt,trials,earlyUpMUAmat)
c = colorbar;
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
ylabel('Ripple event #')
caxis([20 100])
colormap hot
colorbar('off')
subplot(143)
trials = 1:length(lateUp);
imagesc(tt,trials,lateUpMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 100])
colormap hot
colorbar('off')
subplot(142)
trials = 1:length(betweenUp);
imagesc(tt,trials,betweenUpMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 100])
colormap hot
colorbar('off')
subplot(144)
trials = 1:length(upUp);
imagesc(tt,trials,upUpMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 100])
colormap hot
colorbar('off')

%%
lateUpMUAmatctxMUAmean = mean(lateUpMUAmat,1);
lateUpMUAmatctxMUAerror = mad(lateUpMUAmat,1);
hipplateUpMUAmatmean = mean(hipplateUpMUAmat,1);
hipplateUpMUAmaterror = mad(hipplateUpMUAmat,1);
hipplateUpRipplematmean = mean(hipplateUpRipplemat,1);
hipplateUpRipplematerror = mad(hipplateUpRipplemat,1);

figure
hold on
ctx = plot(tt, lateUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,lateUpMUAmatctxMUAmean + lateUpMUAmatctxMUAerror,... 
    lateUpMUAmatctxMUAmean - lateUpMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hipplateUpMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hipplateUpMUAmatmean + hipplateUpMUAmaterror,... 
    hipplateUpMUAmatmean - hipplateUpMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 70])
%%
earlyUpMUAmatctxMUAmean = mean(earlyUpMUAmat,1);
earlyUpMUAmatctxMUAerror = mad(earlyUpMUAmat,1);
hippearlyUpMUAmatmean = mean(hippearlyUpMUAmat,1);
hippearlyUpMUAmaterror = mad(hippearlyUpMUAmat,1);
hippearlyUpRipplematmean = mean(hippearlyUpRipplemat,1);
hippearlyUpRipplematerror = mad(hippearlyUpRipplemat,1);
figure
hold on
ctx = plot(tt, earlyUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,earlyUpMUAmatctxMUAmean + earlyUpMUAmatctxMUAerror,... 
    earlyUpMUAmatctxMUAmean - earlyUpMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hippearlyUpMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippearlyUpMUAmatmean + hippearlyUpMUAmaterror,... 
    hippearlyUpMUAmatmean - hippearlyUpMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])
%%
betweenUpMUAmatctxMUAmean = mean(betweenUpMUAmat,1);
betweenUpMUAmatctxMUAerror = mad(betweenUpMUAmat,1);
hippbetweenUpMUAmatmean = mean(hippbetweenUpMUAmat,1);
hippbetweenUpMUAmaterror = mad(hippbetweenUpMUAmat,1);
hippbetweenUpRipplematmean = mean(hippbetweenUpRipplemat,1);
hippbetweenUpRipplematerror = mad(hippbetweenUpRipplemat,1);
figure
hold on
ctx = plot(tt, betweenUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,betweenUpMUAmatctxMUAmean + betweenUpMUAmatctxMUAerror,... 
    betweenUpMUAmatctxMUAmean - betweenUpMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hippbetweenUpMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippbetweenUpMUAmatmean + hippbetweenUpMUAmaterror,... 
    hippbetweenUpMUAmatmean - hippbetweenUpMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])

%%
upUpMUAmatctxMUAmean = mean(upUpMUAmat,1);
upUpMUAmatctxMUAerror = mad(upUpMUAmat,1);
hippupUpMUAmatmean = mean(hippupUpMUAmat,1);
hippupUpMUAmaterror = mad(hippupUpMUAmat,1);
hippupUpRipplematmean = mean(hippupUpRipplemat,1);
hippupUpRipplematerror = mad(hippupUpRipplemat,1);
figure
hold on
ctx = plot(tt, upUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,upUpMUAmatctxMUAmean + upUpMUAmatctxMUAerror,... 
    upUpMUAmatctxMUAmean - upUpMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hippupUpMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippupUpMUAmatmean + hippupUpMUAmaterror,... 
    hippupUpMUAmatmean - hippupUpMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])
%% Clasification of ripple occuridng during "Down-States"
earlyDown = zeros(size(pksRippleIdx,1),1);
lateDown = zeros(size(pksRippleIdx,1),1);
betweenDown = zeros(size(pksRippleIdx,1),1);
DownDown = zeros(size(pksRippleIdx,1),1);

for kk = 1: size(pksRippleIdx,1)
   if  ctxCentermeddian(kk) < 24 % if ripple occur during upstate
     if   ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) > 24
         betweenDown(kk) = pksRippleIdx(kk);
     elseif ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) < 24
         earlyDown(kk) = pksRippleIdx(kk);
     elseif ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) > 24
         lateDown(kk) = pksRippleIdx(kk);
     elseif ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) < 24
         DownDown(kk) = pksRippleIdx(kk); 
     end
   end
end
 
earlyDown(earlyDown == 0) = [];
lateDown(lateDown == 0) = [];
betweenDown(betweenDown == 0) = [];
DownDown(DownDown == 0) = [];

%% Matrix calculation for Ripples occuring during Down-states
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
earlyDownMUAmat = zeros(size(earlyDown,1),Veclenght);
hippearlyDownMUAmat = zeros(size(earlyDown,1),Veclenght);
hippearlyDownRipplemat = zeros(size(earlyDown,1),Veclenght);
for kk = 1: size(earlyDown,1)
    earlyDownMUAmat(kk,:) = ctxMUAenvpsm(earlyDown(kk) - eventExt : earlyDown(kk)+ eventExt);
    hippearlyDownMUAmat(kk,:) = CA3spkDensity(earlyDown(kk) - eventExt : earlyDown(kk)+ eventExt);
    hippearlyDownRipplemat(kk,:) = ZhippRippleenv(earlyDown(kk) - eventExt : earlyDown(kk)+ eventExt);
end

lateDownMUAmat= zeros(size(lateDown,1),Veclenght);
hipplateDownMUAmat = zeros(size(lateDown,1),Veclenght);
hipplateDownRipplemat = zeros(size(lateDown,1),Veclenght);
for kk = 1: size(lateDown,1)
    lateDownMUAmat(kk,:) = ctxMUAenvpsm(lateDown(kk) - eventExt : lateDown(kk)+ eventExt);
    hipplateDownMUAmat(kk,:) = CA3spkDensity(lateDown(kk) - eventExt : lateDown(kk)+ eventExt);
    hipplateDownRipplemat(kk,:) = ZhippRippleenv(lateDown(kk) - eventExt : lateDown(kk)+ eventExt);
end

betweenDownMUAmat= zeros(size(betweenDown,1),Veclenght);
hippbetweenDownMUAmat = zeros(size(betweenDown,1),Veclenght);
hippbetweenDownRipplemat = zeros(size(betweenDown,1),Veclenght);
for kk = 1: size(betweenDown,1)
    betweenDownMUAmat(kk,:) = ctxMUAenvpsm(betweenDown(kk) - eventExt : betweenDown(kk)+ eventExt);
    hippbetweenDownMUAmat(kk,:) = CA3spkDensity(betweenDown(kk) - eventExt : betweenDown(kk)+ eventExt);
    hippbetweenDownRipplemat(kk,:) = ZhippRippleenv(betweenDown(kk) - eventExt : betweenDown(kk)+ eventExt);
end

DownDownMUAmat= zeros(size(DownDown,1),Veclenght);
hippDownDownMUAmat = zeros(size(DownDown,1),Veclenght);
hippDownDownRipplemat = zeros(size(DownDown,1),Veclenght);
for kk = 1: size(DownDown,1)
    DownDownMUAmat(kk,:) = ctxMUAenvpsm(DownDown(kk) - eventExt : DownDown(kk)+ eventExt);
    hippDownDownMUAmat(kk,:) = CA3spkDensity(DownDown(kk) - eventExt : DownDown(kk)+ eventExt);
    hippDownDownRipplemat(kk,:) = ZhippRippleenv(DownDown(kk) - eventExt : DownDown(kk)+ eventExt);
end
subplot(141)
trials = 1:length(earlyDown);
imagesc(tt,trials,earlyDownMUAmat)
c = colorbar;
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
ylabel('Ripple event #')
caxis([20 100])
colormap hot
colorbar('off')
subplot(143)
trials = 1:length(lateDown);
imagesc(tt,trials,lateDownMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 100])
colormap hot
colorbar('off')
subplot(142)
trials = 1:length(betweenDown);
imagesc(tt,trials,betweenDownMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 100])
colormap hot
colorbar('off')
subplot(144)
trials = 1:length(DownDown);
imagesc(tt,trials,DownDownMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 100])
colormap hot
colorbar('off')
ylim([1 3])

%%
lateDownMUAmatctxMUAmean = mean(lateDownMUAmat,1);
lateDownMUAmatctxMUAerror = mad(lateDownMUAmat,1);
hipplateDownMUAmatmean = mean(hipplateDownMUAmat,1);
hipplateDownMUAmaterror = mad(hipplateDownMUAmat,1);
hipplateDownRipplematmean = mean(hipplateDownRipplemat,1);
hipplateDownRipplematerror = mad(hipplateDownRipplemat,1);

figure
hold on
ctx = plot(tt, lateDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,lateDownMUAmatctxMUAmean + lateDownMUAmatctxMUAerror,... 
    lateDownMUAmatctxMUAmean - lateDownMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hipplateDownMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hipplateDownMUAmatmean + hipplateDownMUAmaterror,... 
    hipplateDownMUAmatmean - hipplateDownMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])
%%
earlyDownMUAmatctxMUAmean = mean(earlyDownMUAmat,1);
earlyDownMUAmatctxMUAerror = mad(earlyDownMUAmat,1);
hippearlyDownMUAmatmean = mean(hippearlyDownMUAmat,1);
hippearlyDownMUAmaterror = mad(hippearlyDownMUAmat,1);
hippearlyDownRipplematmean = mean(hippearlyDownRipplemat,1);
hippearlyDownRipplematerror = mad(hippearlyDownRipplemat,1);
figure
hold on
ctx = plot(tt, earlyDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,earlyDownMUAmatctxMUAmean + earlyDownMUAmatctxMUAerror,... 
    earlyDownMUAmatctxMUAmean - earlyDownMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hippearlyDownMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippearlyDownMUAmatmean + hippearlyDownMUAmaterror,... 
    hippearlyDownMUAmatmean - hippearlyDownMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])
%%
betweenDownMUAmatctxMUAmean = mean(betweenDownMUAmat,1);
betweenDownMUAmatctxMUAerror = mad(betweenDownMUAmat,1);
hippbetweenDownMUAmatmean = mean(hippbetweenDownMUAmat,1);
hippbetweenDownMUAmaterror = mad(hippbetweenDownMUAmat,1);
hippbetweenDownRipplematmean = mean(hippbetweenDownRipplemat,1);
hippbetweenDownRipplematerror = mad(hippbetweenDownRipplemat,1);
figure
hold on
ctx = plot(tt, betweenDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,betweenDownMUAmatctxMUAmean + betweenDownMUAmatctxMUAerror,... 
    betweenDownMUAmatctxMUAmean - betweenDownMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hippbetweenDownMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippbetweenDownMUAmatmean + hippbetweenDownMUAmaterror,... 
    hippbetweenDownMUAmatmean - hippbetweenDownMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])

%%
DownDownMUAmatctxMUAmean = mean(DownDownMUAmat,1);
DownDownMUAmatctxMUAerror = mad(DownDownMUAmat,1);
hippDownDownMUAmatmean = mean(hippDownDownMUAmat,1);
hippDownDownMUAmaterror = mad(hippDownDownMUAmat,1);
hippDownDownRipplematmean = mean(hippDownDownRipplemat,1);
hippDownDownRipplematerror = mad(hippDownDownRipplemat,1);
figure
hold on
ctx = plot(tt, DownDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,DownDownMUAmatctxMUAmean + DownDownMUAmatctxMUAerror,... 
    DownDownMUAmatctxMUAmean - DownDownMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hippDownDownMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippDownDownMUAmatmean + hippDownDownMUAmaterror,... 
    hippDownDownMUAmatmean - hippDownDownMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])

%% CA3 MUA during Up- and Down-States
hippMUAUpSatemat = [hipplateUpMUAmatmean; hippearlyUpMUAmatmean; ...
    hippbetweenUpMUAmatmean; hippupUpMUAmatmean];
hippMUAUpSatemean = mean(hippMUAUpSatemat,1);
hippMUAUpSateerror = mad(hippMUAUpSatemat,1);

hippMUADownSatemat = [hipplateDownMUAmatmean; hippearlyDownMUAmatmean; ...
   hippbetweenDownMUAmatmean; hippDownDownMUAmatmean];
hippMUADownSatemean = mean(hippMUADownSatemat,1);
hippMUADownSateerror = mad(hippMUADownSatemat,1);

figure
hold on
hipMUAUpState = plot(tt,hippMUAUpSatemean,... 
    'Color',[10/255 97/255 171/255], 'LineWidth',2);
jbfill(tt,hippMUAUpSatemean + hippMUAUpSateerror,... 
    hippMUAUpSatemean - hippMUAUpSateerror,... 
    [10/255 97/255 171/255],[10/255 97/255 171/255],0,0.2);

hipMUADownState = plot(tt, hippMUADownSatemean,...
    'Color',[213/255 29/255 91/255], 'LineWidth',2);
jbfill(tt,hippMUADownSatemean + hippMUADownSateerror,... 
    hippMUADownSatemean - hippMUADownSateerror,... 
    [213/255 29/255 91/255],[213/255 29/255 91/255],0,0.2);
line([0 0],[0,48], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('mean CA3 MUA (spk/sec) ')
legend([hipMUAUpState,hipMUADownState],'CA3 MUA during Up-States', ...
    'CA3 MUA during Down-States')
legend boxoff
axis([-inf inf 0 55])

%% Ripple power Up- and Down-States
hippRippleUpSatemat = [hipplateUpRipplematmean; hippearlyUpRipplematmean; ...
    hippbetweenUpRipplematmean; hippupUpRipplematmean];
hippRippleUpSatemean = mean(hippRippleUpSatemat,1);
hippRippleUpSateerror = mad(hippRippleUpSatemat,1);

hippRippleDownSatemat = [hipplateDownRipplematmean; hippearlyDownRipplematmean; ...
   hippbetweenDownRipplematmean; hippDownDownRipplematmean];
hippRippleDownSatemean = mean(hippRippleDownSatemat,1);
hippRippleDownSateerror = mad(hippRippleDownSatemat,1);

figure
hold on
hipRippleUpState = plot(tt,hippRippleUpSatemean,... 
    'Color',[10/255 97/255 171/255], 'LineWidth',2);
jbfill(tt,hippRippleUpSatemean + hippRippleUpSateerror,... 
    hippRippleUpSatemean - hippRippleUpSateerror,... 
    [10/255 97/255 171/255],[10/255 97/255 171/255],0,0.2);

hipRippleDownState = plot(tt, hippRippleDownSatemean,...
    'Color',[213/255 29/255 91/255], 'LineWidth',2);
jbfill(tt,hippRippleDownSatemean + hippRippleDownSateerror,... 
    hippRippleDownSatemean - hippRippleDownSateerror,... 
    [213/255 29/255 91/255],[213/255 29/255 91/255],0,0.2);
line([0 0],[0,8], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('mean Ripple(125-250Hz) Z-Score')
legend([hipRippleUpState,hipRippleDownState],'CA3 Ripple during Up-States', ...
    'CA3 Ripple during Down-States')
legend boxoff
axis([-inf inf 0 9])

%%
totalRipdetected = size(lateUpMUAmat,1) + size(earlyUpMUAmat,1) ...
    + size(betweenUpMUAmat,1) + size(upUpMUAmat,1) + ...
    size(lateDownMUAmat,1) + size(earlyDownMUAmat,1) ...
    + size(betweenDownMUAmat,1) + size(DownDownMUAmat,1);
rippleUpStatePorcentage = (size(lateUpMUAmat,1) + size(earlyUpMUAmat,1) ...
    + size(betweenUpMUAmat,1) + size(upUpMUAmat,1))/totalRipdetected;
rippleDownStatePorcentage = (size(lateDownMUAmat,1) + size(earlyDownMUAmat,1) ...
    + size(betweenDownMUAmat,1) + size(DownDownMUAmat,1))/totalRipdetected;
%%
X = categorical({'"Up-State"','"Down-State"'});
X = reordercats(X,{'"Up-State"','"Down-State"'});
b = bar(X,[rippleUpStatePorcentage*100, rippleDownStatePorcentage*100]);
b.FaceColor = 'flat';
b.CData(2,:) = [213/255 29/255 91/255];
b.CData(1,:) = [10/255 1/255 171/255];
ylabel('Ripple Occurrence %')
ylim([0 100])
box off
%%
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
ylabel('Ripple Occurrence %')
ylim([0 65])
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
ylabel('Ripple Occurrence %')
ylim([0 65])
box off
%% verify the code
close all
% h= histogram(ZhippRippleenv(pksRippleIdx),'Normalization','probability',... 
%     'DisplayStyle','stairs', 'LineWidth',1);
% h.BinEdges = 0:15;
% h.BinWidth = 0.4;
% h.EdgeColor = 'k';
% xlabel('Z-Score Ripple Power')
% ylabel('Probability Density')
% legend('Ripple Power')
% legend('boxoff')  
% box off

figure
h= histogram(ctxMUA,'Normalization','probability',... 
    'DisplayStyle','stair','EdgeColor','red', 'LineWidth',1);
% h.BinEdges = 0:200;
h.BinWidth = 2;
% h.EdgeColor = 'k';
xlabel('MUA Ctx ennp. (uV)')
ylabel('Probability Density')
legend('Ctx MUA at Ripple Peak')
legend('boxoff')  
box off

%% Visual Inspection of Cortical activity based on ripple detection
close all
eventExt = 0.1 / sint;
for kk = 650:680 %: size(,1)
   figure
   subplot(611)
   plot(time(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt), ...
       CTXLFP(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),'k')
   axis([-inf inf, -inf inf])
   
   subplot(612)
   plot(time(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),...
       ctxMUAenvpsm(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),'k')
   axis([-inf inf, 0 inf])
   
   subplot(613)
   plot(time(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt), ...
       CA3LFP(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),'r')
   axis([-inf inf, -inf inf])
   
   subplot(614)
   plot(time(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt), ...
       hippRipple(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),'r')
   axis([-inf inf, -inf inf])
   
   subplot(615)
   plot(time(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),...
       ZhippRippleenv(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),'r')
   axis([-inf inf, -inf inf])
   
   subplot(616)
    plot(time(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),... 
        CA3spkDensity(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),... 
        'r') 
    xlabel('Time[sec]')
    ylabel('spks/s')
    axis tight
   axis([-inf inf, 0 inf])
end










